import subprocess
import sys
import os

def is_deeptools_installed():
    """Verifica se deeptools è installato nell'ambiente corrente."""
    try:
        subprocess.check_call([sys.executable, "-m", "deeptools", "--version"])
        return True
    except subprocess.CalledProcessError:
        return False

def activate_conda_env(env_name):
    """Attiva l'ambiente conda specificato."""
    activate_script = 'activate' if os.name == 'nt' else 'source activate'
    subprocess.call(f"{activate_script} {env_name}", shell=True)

def run_compute_matrix(reference_point, before_range, after_range, regions_file, signal_file, output_file, sorted_regions_file):
    """Esegue il comando computeMatrix con i parametri dati."""
    cmd = f"computeMatrix reference-point --referencePoint {reference_point} -b {before_range} -a {after_range} -R {regions_file} -S {signal_file} --skipZeros -o {output_file} --outFileSortedRegions {sorted_regions_file}"
    subprocess.call(cmd, shell=True)

def run_plot_heatmap(matrix_file, out_file_name):
    """Esegue il comando plotHeatmap per generare un heatmap."""
    cmd = f"plotHeatmap -m {matrix_file} -out {out_file_name}"
    subprocess.call(cmd, shell=True)

def main():
    # Verifica se deeptools è installato
    if not is_deeptools_installed():
        # Sostituisci 'your_env_name' con il nome del tuo ambiente Conda che contiene deeptools
        activate_conda_env('your_env_name')

    # Scegli tra computeMatrix e plotHeatmap
    choice = input("Vuoi eseguire computeMatrix o plotHeatmap? (cm/ph): ")

    if choice.lower() == 'cm':
        # Richiedi input utente per i parametri
        reference_point = input("Inserisci il reference point (es. TSS): ")
        before_range = int(input("Inserisci il range prima del reference point (es. 3000): "))
        after_range = int(input("Inserisci il range dopo il reference point (es. 3000): "))
        regions_file = input("Inserisci il percorso del file delle regioni (es. path/to/regions.bed): ")
        signal_file = input("Inserisci il percorso del file del segnale (es. path/to/signal.bw): ")
        output_file = input("Inserisci il percorso del file di output (es. path/to/output.gz): ")
        sorted_regions_file = input("Inserisci il percorso del file output delle regioni ordinate  (es. path/to/sorted_regions.bed): ")

        run_compute_matrix(reference_point, before_range, after_range, regions_file, signal_file, output_file, sorted_regions_file)
    
    elif choice.lower() == 'ph':
        # Richiedi input utente per i parametri di plotHeatmap
        matrix_file = input("Inserisci il percorso del file della matrice (es. path/to/matrix.gz): ")
        out_file_name = input("Inserisci il nome del file di output per l'heatmap (es. heatmap.png): ")

        run_plot_heatmap(matrix_file, out_file_name)

if __name__ == "__main__":
    main()
