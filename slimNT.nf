params.outdir = 'output'

// downloads a list of proteome ids and maps them to assembly ids
process getIds{
    publishDir "$params.outdir", mode: 'copy'
    output:
    path '*.txt'

    script:
    """
    wget -qO - 'https://rest.uniprot.org/proteomes/stream?compressed=true&fields=upid%2Corganism%2Corganism_id%2Cgenome_assembly&format=tsv&query=%28*%29' | gzip -d > mapping.txt
    cat <(wget -O - https://proteininformationresource.org/rps/data/current/75/rpg-75.txt) <(wget -O - https://proteininformationresource.org/download/rps/rpg_virus_all/current/rpg-95.txt) > full_proteomes.txt
    grep '^>' full_proteomes.txt | awk '{ sub(/^>/, ""); print \$1 }' > ids.txt


    awk 'NR==FNR{a[\$1];next} \$1 in a' ids.txt mapping.txt > mapped.db
    cp mapped.db mapped.txt
    awk '{print \$NF > "assembly_list.txt"} 1'  mapped.db

    """
}

// downoads a pre-populated list of assembly ids for testing purposes
process testgetDB {
    publishDir "$params.outdir", mode: 'copy'
    output:
    path '*.db'

    script:
    """
    wget -O test.db https://raw.githubusercontent.com/Jgergely11/slimNT/main/test_db.txt 
    """
}

process getGenomes {
    publishDir "$params.outdir/genomes", mode: 'copy', overwrite: true
    input:
    path test_db 
    output: 
    file '*.fna'
    file '*.txt'

    script:
    """
    for i in \$(cat ${test_db}); do
        # Check if the file exists
        if [ ! -f "\$i.zip" ]; then
            # File does not exist, download it
            curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/\$i/download?include_annotation_type=GENOME_FASTA&filename=\$i.zip" -H "Accept: application/zip"
        fi
    done
    
    for i in *.zip; do
        # Extract files, skipping any entry that causes an error
        unzip -p \$i "*.fna" > "\${i%.zip}.fna" || true
    done
    
    # Create a text file listing all .fna files
    ls *.fna > fna_files.txt
    find -name "*.fna" -size 0 > empty_fna_files.txt

    """
}


// process findEmpty {
//     publishDir "$params.outdir/genomes", mode: 'copy', overwrite: true
    
//     input:
//     path query
//     output:
//     file 'empty_fna_files.txt'

//     script:
//     """
//     find -name ${query} size > empty_fna_files.txt
//     """
// }
workflow {
    // Define the workflow
    getIds() // Execute the getIds process
    //testgetDB() // Execute the getDB process
    def files = Channel.fromPath('output/test.db')
    getGenomes(files) // Execute the getGenomes process with the output from getDB
    //def query = Channel.fromPath('/output/genomes/')
    //findEmpty(query) // Execute the findEmpty process in the genomes directory
}
