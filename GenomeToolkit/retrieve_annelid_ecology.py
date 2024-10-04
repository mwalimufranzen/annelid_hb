import requests
import csv

def get_worms_species_info(species_name):
    """
    Retrieve species information from the WoRMS database, including habitat/environment.
    
    Args:
        species_name (str): The genus and species name (e.g., "Aplysia kurodai").

    Returns:
        dict: A dictionary containing species information such as habitat and environment.
    """
    worms_search_url = f"https://www.marinespecies.org/rest/AphiaRecordsByMatchNames?scientificnames[]={species_name}&marine_only=false"
    response = requests.get(worms_search_url)
    species_info = {"species": species_name, "habitat": "N/A", "environment": "N/A"}

    if response.status_code == 200 and response.json():
        # Extract the AphiaID from the first matching record
        records = response.json()
        if records and records[0]:
            aphia_id = records[0][0].get("AphiaID", None)
            if aphia_id:
                # Retrieve detailed species information using the AphiaID
                details_url = f"https://www.marinespecies.org/rest/AphiaRecordByAphiaID/{aphia_id}"
                details_response = requests.get(details_url)
                if details_response.status_code == 200:
                    details_data = details_response.json()
                    # Use the 'environment' field to determine habitat context
                    species_info["Environment"] = details_data.get("environment", "N/A")

    return species_info

def retrieve_annelid_info(species_list_file, output_file):
    """
    Read a species list and retrieve WoRMS data for each species.
    
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
            print(f"Retrieving data for species: {species_name}")
            worms_info = get_worms_species_info(species_name)
            
            # Write the data to the output file
            writer.writerow([species_name, globin_count, avg_length, worms_info["environment"]])

    print(f"Species information written to '{output_file}'.")

# File paths for input species list and output
input_species_file = "species_globin_summary.txt"  # Input file from previous steps
output_species_file = "annelid_ecology_summary_worms.txt"  # Output file to save the ecology information

# Run the script to retrieve information and save it to the output file
retrieve_annelid_info(input_species_file, output_species_file)
