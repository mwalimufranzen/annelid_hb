import requests
import csv

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
    Retrieve species information from the WoRMS database using the API.

    Args:
        species_name (str): The genus and species name (e.g., "Aplysia kurodai").

    Returns:
        dict: A dictionary containing species information such as environment.
    """
    # Search WoRMS API by scientific name to get the AphiaID and environment data
    worms_search_url = f"https://www.marinespecies.org/rest/AphiaRecordsByMatchNames?scientificnames[]={species_name}&marine_only=false"
    response = requests.get(worms_search_url)
    species_info = {"species": species_name, "environment": "N/A"}

    if response.status_code == 200 and response.json():
        records = response.json()
        if records and records[0]:
            record = records[0][0]  # Extract the first match
            species_info["environment"] = get_environment_from_api_response(record)

    return species_info

def retrieve_annelid_info(species_list_file, output_file):
    """
    Read a species list and retrieve WoRMS data for each species using the API.

    Args:
        species_list_file (str): Path to the input species list file.
        output_file (str): Path to the output file for saving species information.
    """
    with open(species_list_file, "r") as infile, open(output_file, "w", newline='') as outfile:
        reader = csv.reader(infile, delimiter="\t")
        writer = csv.writer(outfile, delimiter="\t")

        # Write the header row
        writer.writerow(["Species", "Number_of_Globins", "Average_Length", "Environment"])

        # Skip the header of the input file
        next(reader)

        for row in reader:
            species_name, globin_count, avg_length = row[:3]
            worms_info = get_html_worms_info(species_name)

            # Write the data to the output file
            writer.writerow([species_name, globin_count, avg_length, worms_info["environment"]])

    print(f"Species information written to '{output_file}'.")

# File paths for input species list and output
input_species_file = "species_globin_summary.txt"  # Input file from previous steps
output_species_file = "annelid_ecology_summary_worms_environment.txt"  # Output file to save the ecology information

# Run the script to retrieve information and save it to the output file
retrieve_annelid_info(input_species_file, output_species_file)
