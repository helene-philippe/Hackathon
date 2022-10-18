samples = ['A', 'B']
rule all:
  input:
    expand("stats/{sample}.txt", sample=samples)


list_sra = ["SRR628589"]

rule import_rsa:
  input: list_sra
  output: "~/rsa_files/*.rsa"
  container:
  shell:
    """prefetch $list_sra -O ~/sra_files/"""

rule fastq_dump:
  input: "~/rsa_files/*.rsa"
  output:"~/fastq_files/*.fastq"
  container:
  shell:
    """fastq-dump ~/sra_files/*/*.sra -O ~/fastq_files/""" # comment on tape commandes ? telles quelles ?

rule fastqc:
  input: "~/fastq_files/*.fastq"
  output: ""
  container:
  shell:
    """fastqc ~/fastq_files/*.fastq"""


rule count:
  input: "{sample}.fastq"
  output: "stats/{sample}.txt"
  run:
    count = 0
    with open(input[0], 'r') as fin:
      for line in fin.readlines():
        count += 1
    with open(output[0], "w") as fout:
      N = int(count / 4)
      fout.write("{}".format(N))

rule fastqc:
  input: "{sample}.fastq.gz"
  output:
    html="fastqc/{sample}_fastqc.html",
    zip="fastqc/{sample}_fastqc.zip"
  log: "logs/{sample}.log"
  shell: "fastqc {input} >{log} 2>&1"