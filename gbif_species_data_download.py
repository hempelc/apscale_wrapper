#!/usr/bin/env python3

"""
A script to download GBIF occurrence species data.

Credit: Code modified from TaxonTableTools

By Chris Hempel (chris.hempel@simplexdna.com) on May 23 2025
"""

import pandas as pd
import time
import requests_html
import aiohttp
import asyncio
from aiohttp_retry import RetryClient, ExponentialRetry
from tqdm.asyncio import tqdm


# Function to standardize species names based on GBIF taxonomy
def gbif_parent_check(phylum_name, species_name):
    """
    Standardizes species against the GBIF API (when in doubt based on phylum).

    Returns:
        str: The standardized species name or None if no match.
    """
    time.sleep(0.1)
    with requests_html.HTMLSession() as session:
        request_name = "%20".join(species_name.split(" "))
        response = session.get(
            f"https://api.gbif.org/v1/species/match?verbose=true&name={request_name}&limit=1"
        )
        api_response_json = response.json()

        if (
            "note" in api_response_json
            and "Multiple equal matches" in api_response_json["note"]
        ):
            for match in api_response_json.get("alternatives", []):
                if phylum_name == match.get("phylum", None):
                    return match.get("species", None)
        return api_response_json.get("species", None)


# Wrapper function for the standardization of species names based on GBIF taxonomy
def gbif_check_taxonomy(df):
    taxon_table_df = df[["phylum", "species"]]
    # Define excpetions that are no real taxon names and drop them in the df. Also only keep unique species
    exceptions = [
        "Taxonomy unreliable - multiple matching taxa",
        "Taxonomy unreliable - percentage similarity threshold for rank not met",
        "Taxonomy unreliable - bitscore and alignment length threshold not met",
        "No match in database",
        "Unknown in PR2 database",
        "Unknown in BOLD database",
        "Unknown in SILVA database",
        "Unknown in MIDORI2 database",
        "Taxonomy unreliable - confidence threshold not met",
        "No match in database",
    ]
    taxon_table_df = taxon_table_df.replace(exceptions, None).dropna().drop_duplicates()
    checked_species = []
    # Standardize names
    for _, row in taxon_table_df.iterrows():
        phylum_name = row["phylum"]
        species_name = row["species"]
        if checked_species_name := gbif_parent_check(phylum_name, species_name):
            checked_species.append(checked_species_name)
    # Drop contamination species
    contamination_species = [
        "Sus scrofa",
        "Bos taurus",
        "Homo sapiens",
        "Gallus gallus",
        "Canis lupus",
        "Felis catus",
    ]
    return [taxon for taxon in checked_species if taxon not in contamination_species]


# Functions to download GBID specimen locations asynchronously
## Custom exception for handling HTTP 503 errors
class HTTP503Error(Exception):
    pass


## Custom exception for handling SSL errors
class SSLConnectionError(Exception):
    pass


class ServerDisconnectedError(Exception):
    pass


timeout = aiohttp.ClientTimeout(total=120, connect=20, sock_connect=20, sock_read=40)
retry_options = ExponentialRetry(
    attempts=5,
    exceptions={
        HTTP503Error,
        SSLConnectionError,
        ServerDisconnectedError,
    },  # Retry on 3 known error
)


async def fetch_occurrence(retry_session, taxon_name, country_code):
    request_name = "%20".join(taxon_name.split(" "))
    url = f"https://api.gbif.org/v1/occurrence/search?scientificName={request_name}&country={country_code}"
    try:
        async with retry_session.get(url, timeout=timeout) as response:
            if response.status == 200:
                try:
                    api_response_json = await response.json()
                    return api_response_json.get("count", 0)
                except aiohttp.ContentTypeError:
                    print(f"Unexpected content type at URL: {url}")
                    return 0
            elif response.status == 503:
                raise HTTP503Error(
                    f"Service unavailable for {taxon_name} in {country_code}"
                )
            else:
                print(
                    f"Error fetching data for {taxon_name} in {country_code}: HTTP {response.status}"
                )
                return 0
    except aiohttp.client_exceptions.ServerDisconnectedError as e:
        raise ServerDisconnectedError(
            f"Server connection error for {taxon_name} in {country_code}"
        ) from e
    except aiohttp.client_exceptions.ClientOSError as e:
        raise SSLConnectionError(
            f"SSL connection error for {taxon_name} in {country_code}"
        ) from e


async def fetch_all_occurrences(retry_session, taxon_name, country_codes):
    tasks = [
        fetch_occurrence(retry_session, taxon_name, country_code)
        for country_code in country_codes
    ]
    results = await asyncio.gather(*tasks, return_exceptions=True)

    # Process each result to handle exceptions individually
    final_results = []
    for result, country_code in zip(results, country_codes):
        if isinstance(result, HTTP503Error):
            print(
                f"Service unavailable for {taxon_name} in {country_code}. Setting counts to 0."
            )
            final_results.append(0)  # Default value for unavailable data
        elif isinstance(result, SSLConnectionError):
            print(
                f"SSL connection error for {taxon_name} in {country_code}. Setting counts to 0."
            )
            final_results.append(0)  # Default value for SSL errors
        elif isinstance(result, ServerDisconnectedError):
            print(
                f"Server disconnected for {taxon_name} in {country_code}. Setting counts to 0."
            )
            final_results.append(0)  # Default value for server disconnects
        else:
            final_results.append(result)
    return final_results


async def async_main(gbif_standardized_species, country_codes, occurrence_df):
    async with aiohttp.ClientSession() as session:
        async with RetryClient(session, retry_options=retry_options) as retry_session:
            for taxon_name in tqdm(
                gbif_standardized_species,
                desc="Downloading GBIF species location data",
            ):
                occurrence_list = await fetch_all_occurrences(
                    retry_session, taxon_name, country_codes
                )
                occurrence_df[taxon_name] = occurrence_list
    return occurrence_df


# Function to get species occurrence data per country
def gbif_species_data_per_country(gbif_standardized_species_list):

    # Return empty dictionary and None for the plots if the species list is empty, effectively skipping this step
    if not gbif_standardized_species_list:
        return {}, None

    # Define a dictionary with all countries and codes on Earth
    country_codes_dict = {
        "Andorra": ["AD", "Europe", "Palearctic"],
        "United Arab Emirates": ["AE", "Asia", "Palearctic"],
        "Afghanistan": ["AF", "Asia", "Palearctic"],
        "Antigua and Barbuda": ["AG", "North America", "Neotropical"],
        "Anguilla": ["AI", "North America", "Neotropical"],
        "Albania": ["AL", "Europe", "Palearctic"],
        "Armenia": ["AM", "Asia", "Palearctic"],
        "Angola": ["AO", "Africa", "Afrotropical"],
        "Antarctica": ["AQ", "Antarctica", "Antarctic"],
        "Argentina": ["AR", "South America", "Neotropical"],
        "American Samoa": ["AS", "Oceania", "Oceanian"],
        "Austria": ["AT", "Europe", "Palearctic"],
        "Australia": ["AU", "Oceania", "Australasian"],
        "Aruba": ["AW", "North America", "Neotropical"],
        "Åland Islands": ["AX", "Europe", "Palearctic"],
        "Azerbaijan": ["AZ", "Asia", "Palearctic"],
        "Bosnia and Herzegovina": ["BA", "Europe", "Palearctic"],
        "Barbados": ["BB", "North America", "Neotropical"],
        "Bangladesh": ["BD", "Asia", "Indomalayan"],
        "Belgium": ["BE", "Europe", "Palearctic"],
        "Burkina Faso": ["BF", "Africa", "Afrotropical"],
        "Bulgaria": ["BG", "Europe", "Palearctic"],
        "Bahrain": ["BH", "Asia", "Palearctic"],
        "Burundi": ["BI", "Africa", "Afrotropical"],
        "Benin": ["BJ", "Africa", "Afrotropical"],
        "Saint Barthélemy": ["BL", "North America", "Neotropical"],
        "Bermuda": ["BM", "North America", "Neotropical"],
        "Brunei Darussalam": ["BN", "Asia", "Indomalayan"],
        "Bolivia": ["BO", "South America", "Neotropical"],
        "Bonaire, Sint Eustatius and Saba": [
            "BQ",
            "North America",
            "Neotropical",
        ],
        "Brazil": ["BR", "South America", "Neotropical"],
        "Bahamas": ["BS", "North America", "Neotropical"],
        "Bhutan": ["BT", "Asia", "Indomalayan"],
        "Bouvet Island": ["BV", "Antarctica", "Antarctic"],
        "Botswana": ["BW", "Africa", "Afrotropical"],
        "Belarus": ["BY", "Europe", "Palearctic"],
        "Belize": ["BZ", "North America", "Neotropical"],
        "Canada": ["CA", "North America", "Nearctic"],
        "Cocos (Keeling) Islands": ["CC", "Asia", "Indomalayan"],
        "Congo (Democratic Republic)": ["CD", "Africa", "Afrotropical"],
        "Central African Republic": ["CF", "Africa", "Afrotropical"],
        "Congo": ["CG", "Africa", "Afrotropical"],
        "Switzerland": ["CH", "Europe", "Palearctic"],
        "Côte d'Ivoire": ["CI", "Africa", "Afrotropical"],
        "Cook Islands": ["CK", "Oceania", "Oceanian"],
        "Chile": ["CL", "South America", "Neotropical"],
        "Cameroon": ["CM", "Africa", "Afrotropical"],
        "China": ["CN", "Asia", "Palearctic/Indomalayan"],
        "Colombia": ["CO", "South America", "Neotropical"],
        "Costa Rica": ["CR", "North America", "Neotropical"],
        "Cuba": ["CU", "North America", "Neotropical"],
        "Cabo Verde": ["CV", "Africa", "Afrotropical"],
        "Curaçao": ["CW", "North America", "Neotropical"],
        "Christmas Island": ["CX", "Asia", "Indomalayan"],
        "Cyprus": ["CY", "Asia", "Palearctic"],
        "Czechia": ["CZ", "Europe", "Palearctic"],
        "Germany": ["DE", "Europe", "Palearctic"],
        "Djibouti": ["DJ", "Africa", "Afrotropical"],
        "Denmark": ["DK", "Europe", "Palearctic"],
        "Dominica": ["DM", "North America", "Neotropical"],
        "Dominican Republic": ["DO", "North America", "Neotropical"],
        "Algeria": ["DZ", "Africa", "Palearctic"],
        "Ecuador": ["EC", "South America", "Neotropical"],
        "Estonia": ["EE", "Europe", "Palearctic"],
        "Egypt": ["EG", "Africa", "Palearctic"],
        "Western Sahara": ["EH", "Africa", "Palearctic"],
        "Eritrea": ["ER", "Africa", "Afrotropical"],
        "Spain": ["ES", "Europe", "Palearctic"],
        "Ethiopia": ["ET", "Africa", "Afrotropical"],
        "Finland": ["FI", "Europe", "Palearctic"],
        "Fiji": ["FJ", "Oceania", "Australasian"],
        "Falkland Islands": ["FK", "South America", "Neotropical"],
        "Micronesia": ["FM", "Oceania", "Oceanian"],
        "Faroe Islands": ["FO", "Europe", "Palearctic"],
        "France": ["FR", "Europe", "Palearctic"],
        "Gabon": ["GA", "Africa", "Afrotropical"],
        "United Kingdom": ["GB", "Europe", "Palearctic"],
        "Grenada": ["GD", "North America", "Neotropical"],
        "Georgia": ["GE", "Asia", "Palearctic"],
        "French Guiana": ["GF", "South America", "Neotropical"],
        "Guernsey": ["GG", "Europe", "Palearctic"],
        "Ghana": ["GH", "Africa", "Afrotropical"],
        "Gibraltar": ["GI", "Europe", "Palearctic"],
        "Greenland": ["GL", "North America", "Nearctic"],
        "Gambia": ["GM", "Africa", "Afrotropical"],
        "Guinea": ["GN", "Africa", "Afrotropical"],
        "Guadeloupe": ["GP", "North America", "Neotropical"],
        "Equatorial Guinea": ["GQ", "Africa", "Afrotropical"],
        "Greece": ["GR", "Europe", "Palearctic"],
        "South Georgia and the South Sandwich Islands": [
            "GS",
            "Antarctica",
            "Antarctic",
        ],
        "Guatemala": ["GT", "North America", "Neotropical"],
        "Guam": ["GU", "Oceania", "Oceanian"],
        "Guinea-Bissau": ["GW", "Africa", "Afrotropical"],
        "Guyana": ["GY", "South America", "Neotropical"],
        "Hong Kong": ["HK", "Asia", "Indomalayan"],
        "Heard Island and McDonald Islands": ["HM", "Antarctica", "Antarctic"],
        "Honduras": ["HN", "North America", "Neotropical"],
        "Croatia": ["HR", "Europe", "Palearctic"],
        "Haiti": ["HT", "North America", "Neotropical"],
        "Hungary": ["HU", "Europe", "Palearctic"],
        "Indonesia": ["ID", "Asia", "Indomalayan"],
        "Ireland": ["IE", "Europe", "Palearctic"],
        "Israel": ["IL", "Asia", "Palearctic"],
        "Isle of Man": ["IM", "Europe", "Palearctic"],
        "India": ["IN", "Asia", "Indomalayan"],
        "British Indian Ocean Territory": ["IO", "Asia", "Indomalayan"],
        "Iraq": ["IQ", "Asia", "Palearctic"],
        "Iran": ["IR", "Asia", "Palearctic"],
        "Iceland": ["IS", "Europe", "Palearctic"],
        "Italy": ["IT", "Europe", "Palearctic"],
        "Jersey": ["JE", "Europe", "Palearctic"],
        "Jamaica": ["JM", "North America", "Neotropical"],
        "Jordan": ["JO", "Asia", "Palearctic"],
        "Japan": ["JP", "Asia", "Palearctic"],
        "Kenya": ["KE", "Africa", "Afrotropical"],
        "Kyrgyzstan": ["KG", "Asia", "Palearctic"],
        "Cambodia": ["KH", "Asia", "Indomalayan"],
        "Kiribati": ["KI", "Oceania", "Oceanian"],
        "Comoros": ["KM", "Africa", "Afrotropical"],
        "Saint Kitts and Nevis": ["KN", "North America", "Neotropical"],
        "Korea (Democratic People's Republic)": ["KP", "Asia", "Palearctic"],
        "Korea (Republic)": ["KR", "Asia", "Palearctic"],
        "Kuwait": ["KW", "Asia", "Palearctic"],
        "Cayman Islands": ["KY", "North America", "Neotropical"],
        "Kazakhstan": ["KZ", "Asia", "Palearctic"],
        "Lao People's Democratic Republic": ["LA", "Asia", "Indomalayan"],
        "Lebanon": ["LB", "Asia", "Palearctic"],
        "Saint Lucia": ["LC", "North America", "Neotropical"],
        "Liechtenstein": ["LI", "Europe", "Palearctic"],
        "Sri Lanka": ["LK", "Asia", "Indomalayan"],
        "Liberia": ["LR", "Africa", "Afrotropical"],
        "Lesotho": ["LS", "Africa", "Afrotropical"],
        "Lithuania": ["LT", "Europe", "Palearctic"],
        "Luxembourg": ["LU", "Europe", "Palearctic"],
        "Latvia": ["LV", "Europe", "Palearctic"],
        "Libya": ["LY", "Africa", "Palearctic"],
        "Morocco": ["MA", "Africa", "Palearctic"],
        "Monaco": ["MC", "Europe", "Palearctic"],
        "Moldova (the Republic of)": ["MD", "Europe", "Palearctic"],
        "Montenegro": ["ME", "Europe", "Palearctic"],
        "Saint Martin (French part)": ["MF", "North America", "Neotropical"],
        "Madagascar": ["MG", "Africa", "Afrotropical"],
        "Marshall Islands": ["MH", "Oceania", "Oceanian"],
        "Republic of North Macedonia": ["MK", "Europe", "Palearctic"],
        "Mali": ["ML", "Africa", "Afrotropical"],
        "Myanmar": ["MM", "Asia", "Indomalayan"],
        "Mongolia": ["MN", "Asia", "Palearctic"],
        "Macao": ["MO", "Asia", "Indomalayan"],
        "Northern Mariana Islands": ["MP", "Oceania", "Oceanian"],
        "Martinique": ["MQ", "North America", "Neotropical"],
        "Mauritania": ["MR", "Africa", "Palearctic"],
        "Montserrat": ["MS", "North America", "Neotropical"],
        "Malta": ["MT", "Europe", "Palearctic"],
        "Mauritius": ["MU", "Africa", "Afrotropical"],
        "Maldives": ["MV", "Asia", "Indomalayan"],
        "Malawi": ["MW", "Africa", "Afrotropical"],
        "Mexico": ["MX", "North America", "Nearctic/Neotropical"],
        "Malaysia": ["MY", "Asia", "Indomalayan"],
        "Mozambique": ["MZ", "Africa", "Afrotropical"],
        "Namibia": ["NA", "Africa", "Afrotropical"],
        "New Caledonia": ["NC", "Oceania", "Australasian"],
        "Niger": ["NE", "Africa", "Afrotropical"],
        "Norfolk Island": ["NF", "Oceania", "Australasian"],
        "Nigeria": ["NG", "Africa", "Afrotropical"],
        "Nicaragua": ["NI", "North America", "Neotropical"],
        "Netherlands": ["NL", "Europe", "Palearctic"],
        "Norway": ["NO", "Europe", "Palearctic"],
        "Nepal": ["NP", "Asia", "Indomalayan"],
        "Nauru": ["NR", "Oceania", "Oceanian"],
        "Niue": ["NU", "Oceania", "Oceanian"],
        "New Zealand": ["NZ", "Oceania", "Australasian"],
        "Oman": ["OM", "Asia", "Palearctic"],
        "Panama": ["PA", "North America", "Neotropical"],
        "Peru": ["PE", "South America", "Neotropical"],
        "French Polynesia": ["PF", "Oceania", "Oceanian"],
        "Papua New Guinea": ["PG", "Oceania", "Australasian"],
        "Philippines": ["PH", "Asia", "Indomalayan"],
        "Pakistan": ["PK", "Asia", "Palearctic"],
        "Poland": ["PL", "Europe", "Palearctic"],
        "Saint Pierre and Miquelon": ["PM", "North America", "Nearctic"],
        "Pitcairn": ["PN", "Oceania", "Oceanian"],
        "Puerto Rico": ["PR", "North America", "Neotropical"],
        "Palestine, State of": ["PS", "Asia", "Palearctic"],
        "Portugal": ["PT", "Europe", "Palearctic"],
        "Palau": ["PW", "Oceania", "Oceanian"],
        "Paraguay": ["PY", "South America", "Neotropical"],
        "Qatar": ["QA", "Asia", "Palearctic"],
        "Réunion": ["RE", "Africa", "Afrotropical"],
        "Romania": ["RO", "Europe", "Palearctic"],
        "Serbia": ["RS", "Europe", "Palearctic"],
        "Russian Federation": ["RU", "Europe/Asia", "Palearctic"],
        "Rwanda": ["RW", "Africa", "Afrotropical"],
        "Saudi Arabia": ["SA", "Asia", "Palearctic"],
        "Solomon Islands": ["SB", "Oceania", "Australasian"],
        "Seychelles": ["SC", "Africa", "Afrotropical"],
        "Sudan": ["SD", "Africa", "Afrotropical"],
        "Sweden": ["SE", "Europe", "Palearctic"],
        "Singapore": ["SG", "Asia", "Indomalayan"],
        "Saint Helena, Ascension and Tristan da Cunha": [
            "SH",
            "Africa",
            "Afrotropical",
        ],
        "Slovenia": ["SI", "Europe", "Palearctic"],
        "Svalbard and Jan Mayen": ["SJ", "Europe", "Palearctic"],
        "Slovakia": ["SK", "Europe", "Palearctic"],
        "Sierra Leone": ["SL", "Africa", "Afrotropical"],
        "San Marino": ["SM", "Europe", "Palearctic"],
        "Senegal": ["SN", "Africa", "Afrotropical"],
        "Somalia": ["SO", "Africa", "Afrotropical"],
        "Suriname": ["SR", "South America", "Neotropical"],
        "South Sudan": ["SS", "Africa", "Afrotropical"],
        "Sao Tome and Principe": ["ST", "Africa", "Afrotropical"],
        "El Salvador": ["SV", "North America", "Neotropical"],
        "Syrian Arab Republic": ["SY", "Asia", "Palearctic"],
        "Eswatini": ["SZ", "Africa", "Afrotropical"],
        "Turks and Caicos Islands": ["TC", "North America", "Neotropical"],
        "Chad": ["TD", "Africa", "Afrotropical"],
        "French Southern Territories": ["TF", "Antarctica", "Antarctic"],
        "Togo": ["TG", "Africa", "Afrotropical"],
        "Thailand": ["TH", "Asia", "Indomalayan"],
        "Tajikistan": ["TJ", "Asia", "Palearctic"],
        "Tokelau": ["TK", "Oceania", "Oceanian"],
        "Timor-Leste": ["TL", "Asia", "Indomalayan"],
        "Turkmenistan": ["TM", "Asia", "Palearctic"],
        "Tunisia": ["TN", "Africa", "Palearctic"],
        "Tonga": ["TO", "Oceania", "Oceanian"],
        "Turkey": ["TR", "Asia", "Palearctic"],
        "Trinidad and Tobago": ["TT", "North America", "Neotropical"],
        "Tuvalu": ["TV", "Oceania", "Oceanian"],
        "Taiwan": ["TW", "Asia", "Indomalayan"],
        "Tanzania": ["TZ", "Africa", "Afrotropical"],
        "Ukraine": ["UA", "Europe", "Palearctic"],
        "Uganda": ["UG", "Africa", "Afrotropical"],
        "United States Minor Outlying Islands": ["UM", "Oceania", "Oceanian"],
        "United States of America": ["US", "North America", "Nearctic"],
        "Uruguay": ["UY", "South America", "Neotropical"],
        "Uzbekistan": ["UZ", "Asia", "Palearctic"],
        "Holy See": ["VA", "Europe", "Palearctic"],
        "Saint Vincent and the Grenadines": [
            "VC",
            "North America",
            "Neotropical",
        ],
        "Venezuela (Bolivarian Republic of)": [
            "VE",
            "South America",
            "Neotropical",
        ],
        "Virgin Islands (British)": ["VG", "North America", "Neotropical"],
        "Virgin Islands (U.S.)": ["VI", "North America", "Neotropical"],
        "Viet Nam": ["VN", "Asia", "Indomalayan"],
        "Vanuatu": ["VU", "Oceania", "Oceanian"],
        "Wallis and Futuna": ["WF", "Oceania", "Oceanian"],
        "Samoa": ["WS", "Oceania", "Oceanian"],
        "Yemen": ["YE", "Asia", "Palearctic"],
        "Mayotte": ["YT", "Africa", "Afrotropical"],
        "South Africa": ["ZA", "Africa", "Afrotropical"],
        "Zambia": ["ZM", "Africa", "Afrotropical"],
        "Zimbabwe": ["ZW", "Africa", "Afrotropical"],
    }

    # Extract country codes from the dictionary keys
    country_codes = [values[0] for values in country_codes_dict.values()]

    # Make a df template
    occurrence_df = pd.DataFrame(
        [
            {"Country": country, "Continent": values[1], "Realm": values[2]}
            for country, values in country_codes_dict.items()
        ]
    )

    # Run the asynchronous GBIF specimen location retrieval function
    asyncio.run(
        async_main(gbif_standardized_species_list, country_codes, occurrence_df)
    )

    return occurrence_df


def download_gbif_species_data(apscale_result_df):
    print("Standardizing species names based on GBIF...")
    gbif_standardized_species_list = gbif_check_taxonomy(apscale_result_df)
    return gbif_species_data_per_country(gbif_standardized_species_list)
