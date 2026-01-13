#!/usr/bin/env python3
"""
AWS HealthOmics Workshop - Master Script
Run the complete workshop from start to finish
"""

import sys
import time

def print_banner(text):
    """Print a formatted banner"""
    print("\n" + "="*70)
    print(f"  {text}")
    print("="*70 + "\n")

def print_section(number, title):
    """Print a section header"""
    print("\n" + "-"*70)
    print(f"  SECTION {number}: {title}")
    print("-"*70 + "\n")

def run_workshop():
    """Execute the complete HealthOmics workshop"""
    
    print_banner("AWS HEALTHOMICS COMPLETE WORKSHOP")
    print("This script will guide you through all workshop modules")
    print("from basic setup to advanced variant querying.\n")
    
    print("Workshop Modules:")
    print("  1. Create HealthOmics Stores")
    print("  2. Upload Data (Reference & Sequences)")
    print("  3. Workflow Management")
    print("  4. Variant & Annotation Management")
    print("  5. Athena Query Examples")
    
    print("\n" + "="*70 + "\n")
    
    # Module 1: Store Creation
    print_section(1, "CREATE HEALTHOMICS STORES")
    print("Creating Sequence, Reference, Variant, and Annotation stores...")
    print("\nTo run this module:")
    print("  python3 01_create_stores.py")
    print("\nThis creates the foundation for storing all your genomic data.")
    input("\nPress Enter to continue to Module 2...")
    
    # Module 2: Data Upload
    print_section(2, "UPLOAD DATA")
    print("Uploading reference genomes and sequencing data...")
    print("\nBefore running this module:")
    print("  1. Upload your FASTA reference to S3")
    print("  2. Upload your FASTQ files to S3")
    print("  3. Note the S3 URIs")
    print("\nTo run this module:")
    print("  python3 02_upload_data.py")
    print("\nWhat happens:")
    print("  • Reference genome is imported and indexed")
    print("  • FASTQ files are compressed (50-70% savings)")
    print("  • Data is encrypted and stored optimally")
    input("\nPress Enter to continue to Module 3...")
    
    # Module 3: Workflows
    print_section(3, "WORKFLOW MANAGEMENT")
    print("Running bioinformatics analysis pipelines...")
    print("\nCommon workflows:")
    print("  • Alignment: BWA-MEM, Bowtie2")
    print("  • Variant Calling: GATK HaplotypeCaller")
    print("  • RNA-Seq: STAR, Salmon")
    print("  • QC: FastQC, MultiQC")
    print("\nTo run this module:")
    print("  python3 03_workflow_management.py")
    print("\nWorkflow execution:")
    print("  1. Define workflow (Nextflow/WDL)")
    print("  2. Configure parameters")
    print("  3. Start run")
    print("  4. Monitor progress")
    print("  5. Retrieve outputs")
    input("\nPress Enter to continue to Module 4...")
    
    # Module 4: Variants
    print_section(4, "VARIANT & ANNOTATION MANAGEMENT")
    print("Loading and managing variant data...")
    print("\nThis module covers:")
    print("  • Creating Variant Stores")
    print("  • Importing VCF files")
    print("  • Loading annotations (ClinVar, dbSNP)")
    print("  • Variant normalization")
    print("\nTo run this module:")
    print("  python3 04_variant_management.py")
    print("\nBenefits:")
    print("  • Queryable variant database")
    print("  • SQL-based analysis")
    print("  • Join with annotations")
    print("  • Cohort analysis capabilities")
    input("\nPress Enter to continue to Module 5...")
    
    # Module 5: Athena Queries
    print_section(5, "ATHENA QUERIES")
    print("Querying variants with SQL...")
    print("\nQuery examples:")
    print("  • Find pathogenic variants")
    print("  • Filter by population frequency")
    print("  • Gene-specific analysis")
    print("  • Genomic region queries")
    print("\nTo run this module:")
    print("  python3 05_athena_queries.py")
    print("\nSample query:")
    print("""
    SELECT v.contigname, v.start, a.gene, a.clinical_significance
    FROM variant_store v
    JOIN annotation_store a ON v.contigname = a.chr AND v.start = a.pos
    WHERE a.clinical_significance = 'Pathogenic'
    AND a.gene IN ('BRCA1', 'BRCA2', 'TP53');
    """)
    input("\nPress Enter to see completion summary...")
    
    # Completion Summary
    print_banner("WORKSHOP COMPLETE!")
    print("Congratulations! You've completed the AWS HealthOmics workshop.")
    print("\nWhat you've learned:")
    print("  ✓ How to structure genomic data in HealthOmics")
    print("  ✓ How to upload and manage sequencing data")
    print("  ✓ How to run bioinformatics workflows at scale")
    print("  ✓ How to load and query variant data")
    print("  ✓ How to perform SQL-based variant analysis")
    
    print("\nNext steps:")
    print("  1. Apply these concepts to your own data")
    print("  2. Customize workflows for your use cases")
    print("  3. Automate pipelines with Step Functions")
    print("  4. Scale to cohort and population studies")
    
    print("\nResources:")
    print("  • Documentation: https://docs.aws.amazon.com/omics/")
    print("  • Sample workflows: https://github.com/aws-samples/amazon-omics-tutorials")
    print("  • Best practices: Check the workshop_guide.md file")
    
    print("\n" + "="*70)
    print("Thank you for completing the AWS HealthOmics Workshop!")
    print("="*70 + "\n")


def main():
    """Main entry point"""
    try:
        run_workshop()
    except KeyboardInterrupt:
        print("\n\nWorkshop interrupted. You can resume at any time.")
        print("Each module can be run independently.")
        sys.exit(0)
    except Exception as e:
        print(f"\nError: {e}")
        print("Please check the troubleshooting section in README.md")
        sys.exit(1)


if __name__ == '__main__':
    main()