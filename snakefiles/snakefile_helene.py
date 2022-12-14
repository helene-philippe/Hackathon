
LIST_SRA = ['SRR628589']

# comment on fait pour créer les répertoires avant de faire les rules ? 
mkdir rsa_files
mkdir fastq_files
mkdir fastqc_files
mkdir logs

# commencer se connecter à VM ?
# On pull image dans chaque rule ?

# quand est-ce qu'on doit préciser expand() dans les output/input des autres rules ?


rule all:
  input: expand("fastqc_files/{list_sra}_fastqc.html", list_sra=LIST_SRA) 
  # expand permet de gérer le fait qu'on veut plusieurs fichiers à la fin
  # le list_sra=LIST_SRA défini ici se réutilise partout après ?


rule import_rsa:
  input: {list_sra}
  output:
      sra_file = "~/rsa_files/{list_sra}.sra"
  container: "docker://pegi3s/sratoolkit:2.10.0"
  shell: "prefetch {input} -O ~/sra_files/"

rule fastq_dump:
  input: "~/rsa_files/{list_sra}.sra"
  output:
      fastq_file="~/fastq_files/{list_sra}.fastq"
  container: "docker://pegi3s/sratoolkit:2.10.0"
  shell:
    """fastq-dump ~/sra_files/{list_sra}/{list_sra}.sra -O ~/fastq_files/""" # comment on tape commandes ? telles quelles ?

# c'est quoi histoires de --nogroup?
# threads sert à quoi ?
rule fastqc:
  input: "~/fastq_files/{list_sra}.fastq"
  output:
      html="~/fastqs_files/{list_sra}_fastqc.html"
  log: ”logs/{list_sra}.log” # on met les logs dans un fichier au lieu de les afficher sinon ça pollue
  container: "docker://biocontainers/fastqc:v0.11.9_cv8"
  shell:
    """
    fastqc ~/fastq_files/{list_sra}.fastq -O ~/fastqc_files/
    sudo rm ~/fastqc_files/{list_sra}.zip
    """
