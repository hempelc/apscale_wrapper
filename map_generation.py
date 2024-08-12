#!/usr/bin/env python3

"""
A script to generate GBIF-based species occurence maps from an ESV or OTU table.

Credit: Code modified from TaxonTableTools

By Chris Hempel (chris.hempel@simplexdna.com) on Jul 17 2024
"""

import pandas as pd
import time
import requests_html
import plotly.express as px
import sys
import os
import aiohttp
import asyncio
import plotly.colors as pc
from tqdm.asyncio import tqdm


# Define a function to standardize species names based on GBIF taxonomy
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


# Define a wrapper function for the standardization of species names based on GBIF taxonomy
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
    # Dropcontamination species
    contamination_species = [
        "Sus scrofa",
        "Bos taurus",
        "Homo sapiens",
        "Gallus gallus",
        "Canis lupus",
        "Felis catus",
    ]
    return [taxon for taxon in checked_species if taxon not in contamination_species]


# Define functions to download GBID specimen locations asynchronously
async def fetch_occurrences(session, taxon_name, country_code):
    request_name = "%20".join(taxon_name.split(" "))
    url = f"https://api.gbif.org/v1/occurrence/search?scientificName={request_name}&country={country_code}"
    async with session.get(url) as response:
        api_response_json = await response.json()
        return api_response_json.get("count", 0)


async def fetch_all_occurrences(session, taxon_name, country_codes):
    tasks = [
        fetch_occurrences(session, taxon_name, country_code)
        for country_code in country_codes
    ]
    return await asyncio.gather(*tasks)


async def async_main(gbif_standardized_species, country_codes, occurrence_df):
    async with aiohttp.ClientSession() as session:
        for taxon_name in tqdm(
            gbif_standardized_species, desc="Downloading GBIF species location data"
        ):
            occurrence_list = await fetch_all_occurrences(
                session, taxon_name, country_codes
            )
            occurrence_df[taxon_name] = occurrence_list
    return occurrence_df


def map_generation(gbif_standardized_species_list):
    """
    Returns:
        species_maps, continent_occurrence_plot
    """
    # Define a dictionary with all countries and codes on Earth
    country_codes_dict = {
        "Andorra": ["AD", "Europe"],
        "United Arab Emirates": ["AE", "Asia"],
        "Afghanistan": ["AF", "Asia"],
        "Antigua and Barbuda": ["AG", "North America"],
        "Anguilla": ["AI", "North America"],
        "Albania": ["AL", "Europe"],
        "Armenia": ["AM", "Asia"],
        "Angola": ["AO", "Africa"],
        "Antarctica": ["AQ", "Antarctica"],
        "Argentina": ["AR", "South America"],
        "American Samoa": ["AS", "Oceania"],
        "Austria": ["AT", "Europe"],
        "Australia": ["AU", "Oceania"],
        "Aruba": ["AW", "North America"],
        "Åland Islands": ["AX", "Europe"],
        "Azerbaijan": ["AZ", "Asia"],
        "Bosnia and Herzegovina": ["BA", "Europe"],
        "Barbados": ["BB", "North America"],
        "Bangladesh": ["BD", "Asia"],
        "Belgium": ["BE", "Europe"],
        "Burkina Faso": ["BF", "Africa"],
        "Bulgaria": ["BG", "Europe"],
        "Bahrain": ["BH", "Asia"],
        "Burundi": ["BI", "Africa"],
        "Benin": ["BJ", "Africa"],
        "Saint Barthélemy": ["BL", "North America"],
        "Bermuda": ["BM", "North America"],
        "Brunei Darussalam": ["BN", "Asia"],
        "Bolivia": ["BO", "South America"],
        "Bonaire, Sint Eustatius and Saba": ["BQ", "North America"],
        "Brazil": ["BR", "South America"],
        "Bahamas": ["BS", "North America"],
        "Bhutan": ["BT", "Asia"],
        "Bouvet Island": ["BV", "Antarctica"],
        "Botswana": ["BW", "Africa"],
        "Belarus": ["BY", "Europe"],
        "Belize": ["BZ", "North America"],
        "Canada": ["CA", "North America"],
        "Cocos (Keeling) Islands": ["CC", "Asia"],
        "Congo (Democratic Republic)": ["CD", "Africa"],
        "Central African Republic": ["CF", "Africa"],
        "Congo": ["CG", "Africa"],
        "Switzerland": ["CH", "Europe"],
        "Côte d'Ivoire": ["CI", "Africa"],
        "Cook Islands": ["CK", "Oceania"],
        "Chile": ["CL", "South America"],
        "Cameroon": ["CM", "Africa"],
        "China": ["CN", "Asia"],
        "Colombia": ["CO", "South America"],
        "Costa Rica": ["CR", "North America"],
        "Cuba": ["CU", "North America"],
        "Cabo Verde": ["CV", "Africa"],
        "Curaçao": ["CW", "North America"],
        "Christmas Island": ["CX", "Asia"],
        "Cyprus": ["CY", "Asia"],
        "Czechia": ["CZ", "Europe"],
        "Germany": ["DE", "Europe"],
        "Djibouti": ["DJ", "Africa"],
        "Denmark": ["DK", "Europe"],
        "Dominica": ["DM", "North America"],
        "Dominican Republic": ["DO", "North America"],
        "Algeria": ["DZ", "Africa"],
        "Ecuador": ["EC", "South America"],
        "Estonia": ["EE", "Europe"],
        "Egypt": ["EG", "Africa"],
        "Western Sahara": ["EH", "Africa"],
        "Eritrea": ["ER", "Africa"],
        "Spain": ["ES", "Europe"],
        "Ethiopia": ["ET", "Africa"],
        "Finland": ["FI", "Europe"],
        "Fiji": ["FJ", "Oceania"],
        "Falkland Islands": ["FK", "South America"],
        "Micronesia": ["FM", "Oceania"],
        "Faroe Islands": ["FO", "Europe"],
        "France": ["FR", "Europe"],
        "Gabon": ["GA", "Africa"],
        "United Kingdom": ["GB", "Europe"],
        "Grenada": ["GD", "North America"],
        "Georgia": ["GE", "Asia"],
        "French Guiana": ["GF", "South America"],
        "Guernsey": ["GG", "Europe"],
        "Ghana": ["GH", "Africa"],
        "Gibraltar": ["GI", "Europe"],
        "Greenland": ["GL", "North America"],
        "Gambia": ["GM", "Africa"],
        "Guinea": ["GN", "Africa"],
        "Guadeloupe": ["GP", "North America"],
        "Equatorial Guinea": ["GQ", "Africa"],
        "Greece": ["GR", "Europe"],
        "South Georgia and the South Sandwich Islands": ["GS", "Antarctica"],
        "Guatemala": ["GT", "North America"],
        "Guam": ["GU", "Oceania"],
        "Guinea-Bissau": ["GW", "Africa"],
        "Guyana": ["GY", "South America"],
        "Hong Kong": ["HK", "Asia"],
        "Heard Island and McDonald Islands": ["HM", "Antarctica"],
        "Honduras": ["HN", "North America"],
        "Croatia": ["HR", "Europe"],
        "Haiti": ["HT", "North America"],
        "Hungary": ["HU", "Europe"],
        "Indonesia": ["ID", "Asia"],
        "Ireland": ["IE", "Europe"],
        "Israel": ["IL", "Asia"],
        "Isle of Man": ["IM", "Europe"],
        "India": ["IN", "Asia"],
        "British Indian Ocean Territory": ["IO", "Asia"],
        "Iraq": ["IQ", "Asia"],
        "Iran": ["IR", "Asia"],
        "Iceland": ["IS", "Europe"],
        "Italy": ["IT", "Europe"],
        "Jersey": ["JE", "Europe"],
        "Jamaica": ["JM", "North America"],
        "Jordan": ["JO", "Asia"],
        "Japan": ["JP", "Asia"],
        "Kenya": ["KE", "Africa"],
        "Kyrgyzstan": ["KG", "Asia"],
        "Cambodia": ["KH", "Asia"],
        "Kiribati": ["KI", "Oceania"],
        "Comoros": ["KM", "Africa"],
        "Saint Kitts and Nevis": ["KN", "North America"],
        "Korea (Democratic People's Republic)": ["KP", "Asia"],
        "Korea (Republic)": ["KR", "Asia"],
        "Kuwait": ["KW", "Asia"],
        "Cayman Islands": ["KY", "North America"],
        "Kazakhstan": ["KZ", "Asia"],
        "Lao People's Democratic Republic": ["LA", "Asia"],
        "Lebanon": ["LB", "Asia"],
        "Saint Lucia": ["LC", "North America"],
        "Liechtenstein": ["LI", "Europe"],
        "Sri Lanka": ["LK", "Asia"],
        "Liberia": ["LR", "Africa"],
        "Lesotho": ["LS", "Africa"],
        "Lithuania": ["LT", "Europe"],
        "Luxembourg": ["LU", "Europe"],
        "Latvia": ["LV", "Europe"],
        "Libya": ["LY", "Africa"],
        "Morocco": ["MA", "Africa"],
        "Monaco": ["MC", "Europe"],
        "Moldova (the Republic of)": ["MD", "Europe"],
        "Montenegro": ["ME", "Europe"],
        "Saint Martin (French part)": ["MF", "North America"],
        "Madagascar": ["MG", "Africa"],
        "Marshall Islands": ["MH", "Oceania"],
        "Republic of North Macedonia": ["MK", "Europe"],
        "Mali": ["ML", "Africa"],
        "Myanmar": ["MM", "Asia"],
        "Mongolia": ["MN", "Asia"],
        "Macao": ["MO", "Asia"],
        "Northern Mariana Islands": ["MP", "Oceania"],
        "Martinique": ["MQ", "North America"],
        "Mauritania": ["MR", "Africa"],
        "Montserrat": ["MS", "North America"],
        "Malta": ["MT", "Europe"],
        "Mauritius": ["MU", "Africa"],
        "Maldives": ["MV", "Asia"],
        "Malawi": ["MW", "Africa"],
        "Mexico": ["MX", "North America"],
        "Malaysia": ["MY", "Asia"],
        "Mozambique": ["MZ", "Africa"],
        "Namibia": ["NA", "Africa"],
        "New Caledonia": ["NC", "Oceania"],
        "Niger": ["NE", "Africa"],
        "Norfolk Island": ["NF", "Oceania"],
        "Nigeria": ["NG", "Africa"],
        "Nicaragua": ["NI", "North America"],
        "Netherlands": ["NL", "Europe"],
        "Norway": ["NO", "Europe"],
        "Nepal": ["NP", "Asia"],
        "Nauru": ["NR", "Oceania"],
        "Niue": ["NU", "Oceania"],
        "New Zealand": ["NZ", "Oceania"],
        "Oman": ["OM", "Asia"],
        "Panama": ["PA", "North America"],
        "Peru": ["PE", "South America"],
        "French Polynesia": ["PF", "Oceania"],
        "Papua New Guinea": ["PG", "Oceania"],
        "Philippines": ["PH", "Asia"],
        "Pakistan": ["PK", "Asia"],
        "Poland": ["PL", "Europe"],
        "Saint Pierre and Miquelon": ["PM", "North America"],
        "Pitcairn": ["PN", "Oceania"],
        "Puerto Rico": ["PR", "North America"],
        "Palestine, State of": ["PS", "Asia"],
        "Portugal": ["PT", "Europe"],
        "Palau": ["PW", "Oceania"],
        "Paraguay": ["PY", "South America"],
        "Qatar": ["QA", "Asia"],
        "Réunion": ["RE", "Africa"],
        "Romania": ["RO", "Europe"],
        "Serbia": ["RS", "Europe"],
        "Russian Federation": ["RU", "Europe"],
        "Rwanda": ["RW", "Africa"],
        "Saudi Arabia": ["SA", "Asia"],
        "Solomon Islands": ["SB", "Oceania"],
        "Seychelles": ["SC", "Africa"],
        "Sudan": ["SD", "Africa"],
        "Sweden": ["SE", "Europe"],
        "Singapore": ["SG", "Asia"],
        "Saint Helena, Ascension and Tristan da Cunha": ["SH", "Africa"],
        "Slovenia": ["SI", "Europe"],
        "Svalbard and Jan Mayen": ["SJ", "Europe"],
        "Slovakia": ["SK", "Europe"],
        "Sierra Leone": ["SL", "Africa"],
        "San Marino": ["SM", "Europe"],
        "Senegal": ["SN", "Africa"],
        "Somalia": ["SO", "Africa"],
        "Suriname": ["SR", "South America"],
        "South Sudan": ["SS", "Africa"],
        "Sao Tome and Principe": ["ST", "Africa"],
        "El Salvador": ["SV", "North America"],
        "Syrian Arab Republic": ["SY", "Asia"],
        "Eswatini": ["SZ", "Africa"],
        "Turks and Caicos Islands": ["TC", "North America"],
        "Chad": ["TD", "Africa"],
        "French Southern Territories": ["TF", "Antarctica"],
        "Togo": ["TG", "Africa"],
        "Thailand": ["TH", "Asia"],
        "Tajikistan": ["TJ", "Asia"],
        "Tokelau": ["TK", "Oceania"],
        "Timor-Leste": ["TL", "Asia"],
        "Turkmenistan": ["TM", "Asia"],
        "Tunisia": ["TN", "Africa"],
        "Tonga": ["TO", "Oceania"],
        "Turkey": ["TR", "Europe"],
        "Trinidad and Tobago": ["TT", "North America"],
        "Tuvalu": ["TV", "Oceania"],
        "Taiwan": ["TW", "Asia"],
        "Tanzania": ["TZ", "Africa"],
        "Ukraine": ["UA", "Europe"],
        "Uganda": ["UG", "Africa"],
        "United States Minor Outlying Islands": ["UM", "Oceania"],
        "United States of America": ["US", "North America"],
        "Uruguay": ["UY", "South America"],
        "Uzbekistan": ["UZ", "Asia"],
        "Holy See": ["VA", "Europe"],
        "Saint Vincent and the Grenadines": ["VC", "North America"],
        "Venezuela (Bolivarian Republic of)": ["VE", "South America"],
        "Virgin Islands (British)": ["VG", "North America"],
        "Virgin Islands (U.S.)": ["VI", "North America"],
        "Viet Nam": ["VN", "Asia"],
        "Vanuatu": ["VU", "Oceania"],
        "Wallis and Futuna": ["WF", "Oceania"],
        "Samoa": ["WS", "Oceania"],
        "Yemen": ["YE", "Asia"],
        "Mayotte": ["YT", "Africa"],
        "South Africa": ["ZA", "Africa"],
        "Zambia": ["ZM", "Africa"],
        "Zimbabwe": ["ZW", "Africa"],
    }

    # Extract country codes from the dictionary keys
    country_codes = [values[0] for values in country_codes_dict.values()]

    # Make a df template
    occurrence_df = pd.DataFrame(list(country_codes_dict), columns=["Country"])

    # Run the asynchronous GBIF specimen location retrieval function
    asyncio.run(
        async_main(gbif_standardized_species_list, country_codes, occurrence_df)
    )

    # Create custom colour scale with grey for abundances of 0
    custom_colors = [[0, "#d3d3d3"]]  # Grey at position 0
    num_viridis_colors = 10
    viridis_colors = pc.sequential.Viridis[:num_viridis_colors]
    # Calculate the interval step for the remaining colors
    interval_step = 1 / (num_viridis_colors - 1)
    # Append the first Viridis color with a position of 0.00000001
    custom_colors.append([0.00000001, viridis_colors[0]])
    # Append the rest of the Viridis colors at even intervals
    for i in range(1, num_viridis_colors):
        interval_position = i * interval_step
        custom_colors.append([interval_position, viridis_colors[i]])
    # Ensure the last color is exactly at position 1
    custom_colors[-1][0] = 1

    # Generate a map per species
    species_maps = {}
    for species in gbif_standardized_species_list:
        # Create the map for specimen counts
        map_gbif_specimen_counts = px.choropleth(
            occurrence_df,
            locations="Country",
            locationmode="country names",
            hover_name="Country",
            hover_data={
                species: True,
                "Country": False,
            },
            color=species,
            scope="world",
            color_continuous_scale=custom_colors,
        )
        map_gbif_specimen_counts.update_layout(
            title_text=f"{species} - GBIF specimen count by country",
            geo=dict(
                showframe=False,
            ),
            coloraxis_colorbar=dict(
                title="GBIF specimen count",
            ),
        )
        map_gbif_specimen_counts.update_traces(
            hovertemplate="<b>%{hovertext}</b><br>GBIF specimen count: %{customdata[0]}"
        )
        species_maps[species] = map_gbif_specimen_counts

    # Generate continent occurrence plot
    ## Generate binary occurrence data per continent and species
    occurrence_df["Continent"] = occurrence_df["Country"].map(
        lambda x: country_codes_dict.get(x, [None, None])[1]
    )
    continent_df = occurrence_df.drop("Country", axis=1).groupby("Continent").sum()
    continent_df[continent_df > 0] = 1
    continent_df = continent_df.reset_index()

    ## Melt the DataFrame to long format
    continent_df_melted = continent_df.melt(
        id_vars=["Continent"], var_name="Species", value_name="Detected"
    )

    ## Generate the bubble plot using Plotly
    continent_occurrence_plot = px.scatter(
        continent_df_melted,
        x="Continent",
        y="Species",
        size="Detected",
        color="Continent",
        hover_name="Continent",
        size_max=10,
        title="Detected species by continent",
        height=30 * len(continent_df_melted["Species"].unique()),
        width=550,
    )
    continent_occurrence_plot.update_xaxes(tickangle=35)

    return species_maps, continent_occurrence_plot


def main():
    if len(sys.argv) != 3:
        print("Usage: map_generation.py esv/otu_table.csv path/to/output_directory/")
        return

    file = sys.argv[1]
    outdir = sys.argv[2]

    df = pd.read_csv(file)

    print("Standardizing species names based on GBIF...")
    gbif_standardized_species = gbif_check_taxonomy(df)
    species_maps, continent_occurrence_plot = map_generation(gbif_standardized_species)

    for species in species_maps:
        species_maps[species].write_html(
            os.path.join(
                outdir,
                f"map_{species}.html",
            )
        )

    continent_occurrence_plot.write_html(
        os.path.join(
            outdir,
            "continent_occurrence_plot.html",
        )
    )


if __name__ == "__main__":
    main()
