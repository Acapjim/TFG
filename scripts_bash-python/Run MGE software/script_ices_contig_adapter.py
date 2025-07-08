import os
import re
import csv
print(" Script iniciado")

def parse_reference_file(ref_path):
    ref_data = {}
    current_key = None
    with open(ref_path, 'r') as ref_file:
        lines = [line.strip() for line in ref_file if line.strip()]
        i = 0
        while i < len(lines):
            line = lines[i]
            if line.startswith("GCA_") or line.startswith("GCF_"):
                current_key = line
                ref_data[current_key] = []
                i += 1
                while i < len(lines) and not (lines[i].startswith("GCA_") or lines[i].startswith("GCF_")):
                    ref_data[current_key].append(lines[i])
                    i += 1
            else:
                i += 1
    return ref_data

def find_matching_key(filename_id, ref_data):
    for key in ref_data.keys():
        if filename_id in key or key in filename_id:
            return key
    return None

def process_tsv_file(tsv_path, ref_data, output_dir):
    base_name = os.path.basename(tsv_path)
    match = re.search(r'(GCA_\d+\.\d+|GCF_\d+\.\d+)', base_name)
    if not match:
        print(f"No GCA/GCF ID found in {base_name}")
        return

    file_id = match.group(1)
    matching_key = find_matching_key(file_id, ref_data)

    if not matching_key:
        print(f"No reference entries found for {file_id} (searching in keys of .txt)")
        return

    entries = ref_data[matching_key]

    new_rows = []
    with open(tsv_path, 'r', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) < 3:
                new_rows.append(row)
                continue
            contig_match = re.match(r'contig_(\d+)', row[2])
            if contig_match:
                index = int(contig_match.group(1)) - 1
                if 0 <= index < len(entries):
                    new_value = entries[index].split()[0]
                    row[2] = new_value
                else:
                    print(f"Index {index} out of range for {base_name}")
            new_rows.append(row)

    output_path = os.path.join(output_dir, base_name.replace("_detected_ME.tsv", "_renamed.tsv"))
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(new_rows)
    print(f"Procesado: {output_path}")

def main():
    carpeta_tsv = "/home/adria_tfg/tfg/data/Mobile_elements/icescreen_results/icescreen_tsv_kept/"             #  Carpeta con los archivos TSV
    archivo_referencia = "/home/adria_tfg/tfg/data/Mobile_elements/accesion_contig_chromosome_list.txt" # Archivo de referencia
    output_dir = "/home/adria_tfg/tfg/data/Mobile_elements/icescreen_results/icescreen_adapted_tsv"            # Carpeta de salida

    os.makedirs(output_dir, exist_ok=True)
    ref_data = parse_reference_file(archivo_referencia)

    for archivo in os.listdir(carpeta_tsv):
        if archivo.endswith("_detected_ME.tsv"):
            full_path = os.path.join(carpeta_tsv, archivo)
            process_tsv_file(full_path, ref_data, output_dir)

if __name__ == "__main__":
    main()

