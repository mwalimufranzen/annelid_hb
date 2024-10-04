import requests
import csv
from bs4 import BeautifulSoup

def parse_html_descriptive_notes(soup):
    """
    Parses the HTML page for various descriptive notes, including Biology and Distribution.

    Args:
        soup (BeautifulSoup): The parsed HTML content of the WoRMS species page.

    Returns:
        str: The combined descriptive notes if found, otherwise "N/A".
    """
    # Retrieve content from different sections: Biology, Distribution, and other notes
    sections = ["Descriptive notes", "Biology", "Distribution", "Additional information"]

    descriptive_info = []

    for section in sections:
        # Look for section labels in the HTML (e.g., "Descriptive notes:", "Biology:")
        section_label = soup.find("strong", string=f"{section}:")
        if section_label:
            # Get the text content of the following sibling element
            section_text = section_label.find_next_sibling("span")
            if section_text:
                descriptive_info.append(f"{section}: {section_text.text.strip()}")

    # Combine all the found descriptive information
    combined_notes = " | ".join(descriptive_info) if descriptive_info else "N/A"
    print(f"DEBUG: Retrieved descriptive notes from HTML: {combined_notes}")  # Debug statement
    return combined_notes

def get_environment_from_api_response(api_response):
    """
    Extract the environment information from the WoRMS API response.

    Args:
        api_response (dict): The parsed JSON response from the WoRMS API.

    Returns:
        str: A string describing the environment type.
    """
    if api_response:
        environment = []
        if api_response.get('isMarine') == 1:
            environment.append("Marine")
        if api_response.get('isBrackish') == 1:
            environment.append("Brackish")
        if api_response.get('isFreshwater') == 1:
            environment.append("Freshwater")
        if api_response.get('isTerrestrial') == 1:
            environment.append("Terrestrial")

        return ", ".join(environment) if environment else "N/A"
    return "N/A"

def get_html_worms_info(species_name):
    """
    Retrieve species information from the WoRMS database HTML page, including environment and descriptive notes.

    Args:
        species_name (str): The genus and species name (e.g., "Aplysia kurodai").

    Returns:
        dict: A dictionary containing species information such as environment and descriptive notes.
    """
    # Search WoRMS API by scientific name to get the AphiaID and environment data
    worms_search_url = f"https://www.marinespecies.org/rest/AphiaRecordsByMatchNames?scientificnames[]={species_name}&marine_only=false"
    print(f"DEBUG: Searching for species: {species_name} using URL: {worms_search_url}")
    response = requests.get(worms_search_url)
    species_info = {"species": species_name, "environment": "N/A", "descriptive_notes": "N/A"}

    if response.status_code == 200 and response.json():
        # Print the full API response to check the content
        records = response.json()
        print(f"DEBUG: Full API Response for {species_name}: {records}")  # Debug statement
        if records and records[0]:
            record = records[0][0]  # Extract the first match
            species_info["environment"] = get_environment_from_api_response(record)
            aphia_id = record.get("AphiaID", None)
            print(f"DEBUG: Found AphiaID: {aphia_id} for species: {species_name}")  # Debug statement

            if aphia_id:
                # Use the AphiaID to access the HTML page for this species
                details_html_url = f"https://www.marinespecies.org/aphia.php?p=taxdetails&id={aphia_id}"
                html_response = requests.get(details_html_url)

                if html_response.status_code == 200:
                    # Print the HTML content for verification
                    print(f"DEBUG: HTML content for {species_name}:\n{html_response.text[:1000]}")  # Limit output to 1000 characters for readability
                    soup = BeautifulSoup(html_response.text, 'html.parser')

                    # Extract descriptive notes using the HTML parser
                    species_info["descriptive_notes"] = parse_html_descriptive_notes(soup)

    return species_info

def retrieve_annelid_info(species_list_file, output_file):
    """
    Read a species list and retrieve WoRMS data for each species using the API and HTML scraping.

    Args:
        species_list_file (str): Path to the input species list file.
        output_file (str): Path to the output file for saving species information.
    """
    with open(species_list_file, "r") as infile, open(output_file, "w", newline='') as outfile:
        reader = csv.reader(infile, delimiter="\t")
        writer = csv.writer(outfile, delimiter="\t")

        # Write the header row
        writer.writerow(["Species", "Number_of_Globins", "Average_Length", "Environment", "Descriptive_Notes"])

        # Skip the header of the input file
        next(reader)

        for row in reader:
            species_name, globin_count, avg_length = row[:3]
            print(f"DEBUG: Retrieving data for species: {species_name}")
            worms_info = get_html_worms_info(species_name)

            # Write the data to the output file
            writer.writerow([species_name, globin_count, avg_length, worms_info["environment"], worms_info["descriptive_notes"]])
            print(f"DEBUG: Written data for {species_name} - Environment: {worms_info['environment']}, Descriptive Notes: {worms_info['descriptive_notes']}")

    print(f"Species information written to '{output_file}'.")

# File paths for input species list and output
input_species_file = "species_globin_summary.txt"  # Input file from previous steps
output_species_file = "annelid_ecology_summary_worms_additional_info.txt"  # Output file to save the ecology information

# Run the script to retrieve information and save it to the output file
retrieve_annelid_info(input_species_file, output_species_file)
