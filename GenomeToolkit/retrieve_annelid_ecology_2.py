import requests
import csv

def get_gbif_species_info(species_name):
    """
    Retrieve species information from the GBIF database.
    
    Args:
        species_name (str): The genus and species name (e.g., "Aplysia kurodai").

    Returns:
        dict: A dictionary containing species information such as habitat, location, and distribution.
    """
    gbif_url = f"https://api.gbif.org/v1/species/match?name={species_name}"
    response = requests.get(gbif_url)
    species_info = {"species": species_name, "habitat": "N/A", "location": "N/A"}

    if response.status_code == 200:
        data = response.json()
        if "usageKey" in data:
            species_key = data["usageKey"]

            # Retrieve details for the species using its GBIF key
            details_url = f"https://api.gbif.org/v1/species/{species_key}"
            details_response = requests.get(details_url)
            if details_response.status_code == 200:
                details_data = details_response.json()
                if "habitats" in details_data:
                    habitat_info = details_data["habitats"]
                    species_info["habitat"] = ", ".join([h["name"] for h in habitat_info])
                if "vernacularNames" in details_data:
                    locations = [v.get("area") for v in details_data["vernacularNames"] if v.get("area")]
                    species_info["location"] = ", ".join(locations) if locations else "N/A"
    return species_info

def retrieve_annelid_info(species_list_file, output_file):
    """
    Read a species list and retrieve GBIF data for each species.
    
    Args:
        species_list_file (str): Path to the input species list file.
        output_file (str): Path to the output file for saving species information.
    """
    with open(species_list_file, "r") as infile, open(output_file, "w", newline='') as outfile:
        reader = csv.reader(infile, delimiter="\t")
        writer = csv.writer(outfile, delimiter="\t")
        
        # Write the header row
        writer.writerow(["Species", "Number_of_Globins", "Average_Length", "Habitat", "Location"])

        # Skip the header of the input file
        next(reader)

        for row in reader:
            species_name, globin_count, avg_length = row[:3]
            print(f"Retrieving data for species: {species_name}")
            gbif_info = get_gbif_species_info(species_name)
            
            # Write the data to the output file
            writer.writerow([species_name, globin_count, avg_length, gbif_info["habitat"], gbif_info["location"]])

    print(f"Species information written to '{output_file}'.")

# File paths for input species list and output
input_species_file = "species_globin_summary.txt"  # Input file from previous steps
output_species_file = "annelid_ecology_summary.txt"  # Output file to save the ecology information

# Run the script to retrieve information and save it to the output file
retrieve_annelid_info(input_species_file, output_species_file)
