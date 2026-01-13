#!/usr/bin/env python3
"""
Module 5: AWS HealthOmics Athena Queries
This script demonstrates how to query variant and annotation data using Athena
"""

import boto3
import json
import time
from botocore.exceptions import ClientError

class HealthOmicsAthenaQuerier:
    """
    Query variants and annotations using Amazon Athena
    
    Supports:
    - Querying variants from Variant Store
    - Querying annotations from Annotation Store
    - Joining variants with annotations
    - Complex filtering and aggregation
    """
    
    def __init__(self, 
                 config_file='/home/claude/healthomics_config.json',
                 output_location='s3://healthomics-workshop-output/athena-results/'):
        """
        Initialize Athena querier
        
        Args:
            config_file (str): Path to configuration file
            output_location (str): S3 location for query results
        """
        with open(config_file, 'r') as f:
            self.config = json.load(f)
        
        self.region = self.config['region']
        self.athena_client = boto3.client('athena', region_name=self.region)
        self.output_location = output_location
        
        print(f"Initialized Athena Querier")
        print(f"Region: {self.region}")
        print(f"Output Location: {output_location}\n")
    
    def execute_query(self, query, wait_for_completion=True):
        """
        Execute an Athena query
        
        Args:
            query (str): SQL query to execute
            wait_for_completion (bool): Wait for query to complete
            
        Returns:
            str: Query execution ID
        """
        print("="*60)
        print("Executing Athena Query")
        print("="*60)
        print(f"Query:\n{query}\n")
        
        try:
            # Start query execution
            response = self.athena_client.start_query_execution(
                QueryString=query,
                ResultConfiguration={
                    'OutputLocation': self.output_location
                }
            )
            
            query_execution_id = response['QueryExecutionId']
            print(f"✓ Query started successfully!")
            print(f"  Execution ID: {query_execution_id}")
            
            if wait_for_completion:
                print("\nWaiting for query to complete...")
                self._wait_for_query_completion(query_execution_id)
                
                # Get results
                results = self.get_query_results(query_execution_id)
                self._display_results(results)
            
            return query_execution_id
            
        except ClientError as e:
            print(f"✗ Error executing query: {e}")
            raise
    
    def _wait_for_query_completion(self, query_execution_id, max_attempts=60):
        """Wait for query to complete"""
        for attempt in range(max_attempts):
            try:
                response = self.athena_client.get_query_execution(
                    QueryExecutionId=query_execution_id
                )
                
                status = response['QueryExecution']['Status']['State']
                
                print(f"  [{attempt+1}/{max_attempts}] Status: {status}", end='\r')
                
                if status == 'SUCCEEDED':
                    print("\n✓ Query completed successfully!")
                    
                    # Show statistics
                    stats = response['QueryExecution']['Statistics']
                    print(f"\nQuery Statistics:")
                    print(f"  Data Scanned: {stats.get('DataScannedInBytes', 0) / 1024 / 1024:.2f} MB")
                    print(f"  Execution Time: {stats.get('EngineExecutionTimeInMillis', 0) / 1000:.2f} seconds")
                    return
                    
                elif status == 'FAILED':
                    print("\n")
                    error = response['QueryExecution']['Status'].get('StateChangeReason', 'Unknown error')
                    raise Exception(f"Query failed: {error}")
                    
                elif status == 'CANCELLED':
                    print("\n")
                    raise Exception("Query was cancelled")
                
                time.sleep(2)
                
            except ClientError as e:
                if attempt == max_attempts - 1:
                    raise
                time.sleep(2)
        
        raise Exception("Query timed out")
    
    def get_query_results(self, query_execution_id, max_results=100):
        """
        Get query results
        
        Args:
            query_execution_id (str): Query execution ID
            max_results (int): Maximum number of results to return
            
        Returns:
            dict: Query results
        """
        try:
            response = self.athena_client.get_query_results(
                QueryExecutionId=query_execution_id,
                MaxResults=max_results
            )
            
            return response
            
        except ClientError as e:
            print(f"Error getting query results: {e}")
            raise
    
    def _display_results(self, results):
        """Display query results in a formatted table"""
        print("\nQuery Results:")
        print("="*60)
        
        result_set = results.get('ResultSet', {})
        rows = result_set.get('Rows', [])
        
        if not rows:
            print("No results found")
            return
        
        # Get column names from first row
        headers = [col.get('VarCharValue', '') for col in rows[0].get('Data', [])]
        
        # Display headers
        header_line = " | ".join(f"{h:15}" for h in headers)
        print(header_line)
        print("-" * len(header_line))
        
        # Display data rows
        for row in rows[1:]:  # Skip header row
            values = [col.get('VarCharValue', 'NULL') for col in row.get('Data', [])]
            data_line = " | ".join(f"{v:15}" for v in values)
            print(data_line)
        
        print("\n")
    
    # Pre-defined queries for common use cases
    
    def query_variants_by_gene(self, gene_name, variant_store_table):
        """
        Query variants in a specific gene
        
        Args:
            gene_name (str): Gene symbol (e.g., 'BRCA1')
            variant_store_table (str): Athena table name for variant store
        """
        query = f"""
        SELECT 
            contigname as chromosome,
            start as position,
            referenceallele as ref,
            alternatealleles as alt,
            names as rsid
        FROM {variant_store_table}
        WHERE gene = '{gene_name}'
        LIMIT 100;
        """
        
        print(f"Querying variants in gene: {gene_name}\n")
        return self.execute_query(query)
    
    def query_pathogenic_variants(self, variant_store_table, annotation_store_table):
        """
        Query pathogenic variants with annotations
        
        Args:
            variant_store_table (str): Athena table for variant store
            annotation_store_table (str): Athena table for annotation store
        """
        query = f"""
        SELECT 
            v.contigname as chromosome,
            v.start as position,
            v.referenceallele as ref,
            v.alternatealleles as alt,
            a.gene,
            a.consequence,
            a.clinical_significance
        FROM {variant_store_table} v
        JOIN {annotation_store_table} a
        ON v.contigname = a.chr
        AND v.start = a.pos
        AND v.referenceallele = a.ref
        AND v.alternatealleles = a.alt
        WHERE a.clinical_significance IN ('Pathogenic', 'Likely_pathogenic')
        LIMIT 100;
        """
        
        print("Querying pathogenic variants with annotations\n")
        return self.execute_query(query)
    
    def query_rare_variants(self, variant_store_table, annotation_store_table, max_af=0.01):
        """
        Query rare variants (low population frequency)
        
        Args:
            variant_store_table (str): Athena table for variant store
            annotation_store_table (str): Athena table for annotation store
            max_af (float): Maximum allele frequency threshold
        """
        query = f"""
        SELECT 
            v.contigname as chromosome,
            v.start as position,
            v.referenceallele as ref,
            v.alternatealleles as alt,
            a.gene,
            a.population_af as allele_frequency,
            a.consequence
        FROM {variant_store_table} v
        JOIN {annotation_store_table} a
        ON v.contigname = a.chr
        AND v.start = a.pos
        WHERE CAST(a.population_af AS DOUBLE) < {max_af}
        ORDER BY CAST(a.population_af AS DOUBLE)
        LIMIT 100;
        """
        
        print(f"Querying rare variants (AF < {max_af})\n")
        return self.execute_query(query)
    
    def query_variant_counts_by_chromosome(self, variant_store_table):
        """
        Count variants per chromosome
        
        Args:
            variant_store_table (str): Athena table for variant store
        """
        query = f"""
        SELECT 
            contigname as chromosome,
            COUNT(*) as variant_count
        FROM {variant_store_table}
        GROUP BY contigname
        ORDER BY 
            CASE 
                WHEN contigname = 'chrX' THEN 23
                WHEN contigname = 'chrY' THEN 24
                WHEN contigname = 'chrM' THEN 25
                ELSE CAST(REPLACE(contigname, 'chr', '') AS INT)
            END;
        """
        
        print("Counting variants by chromosome\n")
        return self.execute_query(query)
    
    def query_missense_variants_with_high_impact(self, 
                                                 variant_store_table,
                                                 annotation_store_table):
        """
        Query missense variants predicted to have high impact
        
        Args:
            variant_store_table (str): Athena table for variant store
            annotation_store_table (str): Athena table for annotation store
        """
        query = f"""
        SELECT 
            v.contigname as chromosome,
            v.start as position,
            v.referenceallele as ref,
            v.alternatealleles as alt,
            a.gene,
            a.consequence,
            a.clinical_significance
        FROM {variant_store_table} v
        JOIN {annotation_store_table} a
        ON v.contigname = a.chr
        AND v.start = a.pos
        WHERE a.consequence LIKE '%missense%'
        AND a.clinical_significance IN ('Pathogenic', 'Likely_pathogenic', 'Uncertain_significance')
        LIMIT 100;
        """
        
        print("Querying high-impact missense variants\n")
        return self.execute_query(query)
    
    def query_variants_in_genomic_region(self, 
                                        variant_store_table,
                                        chromosome,
                                        start_pos,
                                        end_pos):
        """
        Query variants in a specific genomic region
        
        Args:
            variant_store_table (str): Athena table for variant store
            chromosome (str): Chromosome (e.g., 'chr1')
            start_pos (int): Start position
            end_pos (int): End position
        """
        query = f"""
        SELECT 
            contigname as chromosome,
            start as position,
            referenceallele as ref,
            alternatealleles as alt,
            names as rsid,
            qual as quality
        FROM {variant_store_table}
        WHERE contigname = '{chromosome}'
        AND start BETWEEN {start_pos} AND {end_pos}
        ORDER BY start
        LIMIT 100;
        """
        
        print(f"Querying variants in region: {chromosome}:{start_pos}-{end_pos}\n")
        return self.execute_query(query)
    
    def create_variant_summary_statistics(self, variant_store_table):
        """
        Create summary statistics for variants
        
        Args:
            variant_store_table (str): Athena table for variant store
        """
        query = f"""
        SELECT 
            COUNT(*) as total_variants,
            COUNT(DISTINCT contigname) as chromosomes,
            AVG(CAST(qual AS DOUBLE)) as avg_quality,
            MIN(start) as min_position,
            MAX(start) as max_position
        FROM {variant_store_table};
        """
        
        print("Creating variant summary statistics\n")
        return self.execute_query(query)


def demonstrate_common_queries():
    """
    Demonstrate common Athena queries for variant analysis
    """
    print("="*60)
    print("AWS HealthOmics Workshop - Module 5")
    print("Querying Variants with Athena")
    print("="*60)
    print("\nThis module demonstrates:")
    print("  1. Basic variant queries")
    print("  2. Joining variants with annotations")
    print("  3. Filtering by clinical significance")
    print("  4. Population frequency analysis")
    print("  5. Genomic region queries")
    print("\n" + "="*60 + "\n")
    
    # Initialize querier
    querier = HealthOmicsAthenaQuerier()
    
    print("Common Query Examples:")
    print("="*60 + "\n")
    
    # Example query templates
    queries = [
        {
            'name': 'Find Pathogenic Variants',
            'description': 'Query variants classified as pathogenic or likely pathogenic',
            'use_case': 'Clinical variant interpretation'
        },
        {
            'name': 'Rare Variant Discovery',
            'description': 'Find variants with population frequency < 1%',
            'use_case': 'Novel variant discovery, rare disease research'
        },
        {
            'name': 'Gene-Specific Analysis',
            'description': 'Query all variants in a specific gene',
            'use_case': 'Targeted gene analysis, BRCA screening'
        },
        {
            'name': 'Regional Analysis',
            'description': 'Query variants in a genomic region',
            'use_case': 'GWAS studies, structural variant analysis'
        },
        {
            'name': 'Missense Impact',
            'description': 'Find high-impact missense variants',
            'use_case': 'Functional variant analysis'
        },
        {
            'name': 'Variant Statistics',
            'description': 'Summary statistics across entire dataset',
            'use_case': 'QC, data exploration'
        }
    ]
    
    for i, query in enumerate(queries, 1):
        print(f"{i}. {query['name']}")
        print(f"   Description: {query['description']}")
        print(f"   Use Case: {query['use_case']}")
        print()
    
    print("\n" + "="*60)
    print("Sample Query Templates")
    print("="*60 + "\n")
    
    # Show sample SQL queries
    print("Example 1: Find pathogenic variants in BRCA1")
    print("-" * 60)
    sample_query_1 = """
    SELECT 
        v.contigname, v.start, v.referenceallele, v.alternatealleles,
        a.clinical_significance, a.consequence
    FROM variant_store v
    JOIN annotation_store a 
    ON v.contigname = a.chr AND v.start = a.pos
    WHERE a.gene = 'BRCA1' 
    AND a.clinical_significance = 'Pathogenic';
    """
    print(sample_query_1)
    
    print("\nExample 2: Count variants by chromosome")
    print("-" * 60)
    sample_query_2 = """
    SELECT 
        contigname as chromosome,
        COUNT(*) as variant_count
    FROM variant_store
    GROUP BY contigname
    ORDER BY contigname;
    """
    print(sample_query_2)
    
    print("\nExample 3: Find rare, high-impact variants")
    print("-" * 60)
    sample_query_3 = """
    SELECT 
        v.contigname, v.start, a.gene, a.consequence,
        a.population_af, a.clinical_significance
    FROM variant_store v
    JOIN annotation_store a 
    ON v.contigname = a.chr AND v.start = a.pos
    WHERE CAST(a.population_af AS DOUBLE) < 0.01
    AND a.consequence IN ('stop_gained', 'frameshift_variant', 'missense_variant')
    ORDER BY CAST(a.population_af AS DOUBLE);
    """
    print(sample_query_3)
    
    print("\n" + "="*60)
    print("Query Best Practices")
    print("="*60)
    print("""
1. Use partitioning: Query specific chromosomes when possible
2. Limit results: Always use LIMIT for exploratory queries
3. Use indexes: Leverage chromosome and position for filtering
4. Join efficiently: Filter before joining for better performance
5. Monitor costs: Check data scanned in query statistics
6. Use views: Create views for frequently used complex queries
7. Optimize filters: Put most selective filters first
8. Type casting: Cast strings to numbers for numeric comparisons
    """)
    
    print("\n" + "="*60)
    print("Athena Query Module Complete!")
    print("="*60)
    print("\nYou've learned how to:")
    print("  ✓ Write SQL queries for variant data")
    print("  ✓ Join variants with annotations")
    print("  ✓ Filter by clinical significance")
    print("  ✓ Analyze population frequencies")
    print("  ✓ Query specific genomic regions")
    print("  ✓ Generate summary statistics")
    print("\nNext Steps:")
    print("  1. Run these queries on your actual data")
    print("  2. Export results for further analysis")
    print("  3. Integrate with visualization tools")
    print("  4. Automate reporting workflows")
    print("\n" + "="*60)


if __name__ == '__main__':
    demonstrate_common_queries()