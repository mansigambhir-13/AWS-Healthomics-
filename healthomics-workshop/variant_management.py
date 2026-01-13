#!/usr/bin/env python3
"""
Module 4: AWS HealthOmics Variant and Annotation Management
This script demonstrates how to load and manage variant data and annotations
"""

import boto3
import json
import time
from botocore.exceptions import ClientError

class HealthOmicsVariantManager:
    """
    Manages variant data and annotations in HealthOmics
    
    Supports:
    - Loading VCF files into Variant Store
    - Loading annotations into Annotation Store
    - Querying variants
    - Managing variant metadata
    """
    
    def __init__(self, config_file='/home/claude/healthomics_config.json'):
        """Initialize variant manager"""
        with open(config_file, 'r') as f:
            self.config = json.load(f)
        
        self.region = self.config['region']
        self.omics_client = boto3.client('omics', region_name=self.region)
        
        print(f"Initialized Variant Manager")
        print(f"Region: {self.region}")
        print(f"Reference Store: {self.config.get('reference_store_id', 'N/A')}\n")
    
    def create_variant_store_with_reference(self, 
                                           reference_arn,
                                           name="variant-store"):
        """
        Create a Variant Store with a reference genome
        
        Args:
            reference_arn (str): ARN of the reference genome
            name (str): Name for the variant store
            
        Returns:
            str: Variant store ID
        """
        print("="*60)
        print("Creating Variant Store")
        print("="*60)
        print(f"Name: {name}")
        print(f"Reference: {reference_arn}\n")
        
        try:
            response = self.omics_client.create_variant_store(
                name=name,
                description="Workshop variant store for VCF data",
                reference={
                    'referenceArn': reference_arn
                },
                tags={
                    'Project': self.config.get('workshop_name', 'healthomics-workshop')
                }
            )
            
            store_id = response['id']
            print(f"✓ Variant Store created successfully!")
            print(f"  Store ID: {store_id}")
            print(f"  Name: {response['name']}")
            print(f"  Status: {response['status']}")
            
            # Wait for store to be active
            print("\nWaiting for Variant Store to become active...")
            self._wait_for_variant_store(name)
            print("✓ Variant Store is now ACTIVE")
            
            # Update config
            self.config['variant_store_id'] = store_id
            self._save_config()
            
            return store_id
            
        except ClientError as e:
            print(f"✗ Error creating Variant Store: {e}")
            raise
    
    def import_vcf_to_variant_store(self,
                                    s3_vcf_uri,
                                    variant_store_name,
                                    role_arn,
                                    run_left_normalization=True,
                                    annotation_fields=None):
        """
        Import VCF file into Variant Store
        
        Args:
            s3_vcf_uri (str): S3 URI of the VCF file
            variant_store_name (str): Name of the variant store
            role_arn (str): IAM role ARN for import job
            run_left_normalization (bool): Normalize variants
            annotation_fields (dict): Optional annotation field mappings
            
        Returns:
            str: Import job ID
        """
        print("="*60)
        print("Importing VCF to Variant Store")
        print("="*60)
        print(f"Source: {s3_vcf_uri}")
        print(f"Destination: {variant_store_name}")
        print(f"Normalization: {run_left_normalization}\n")
        
        try:
            # Prepare import request
            import_params = {
                'destinationName': variant_store_name,
                'roleArn': role_arn,
                'items': [{
                    'source': s3_vcf_uri
                }],
                'runLeftNormalization': run_left_normalization
            }
            
            if annotation_fields:
                import_params['annotationFields'] = annotation_fields
            
            # Start import job
            response = self.omics_client.start_variant_import_job(**import_params)
            
            job_id = response['id']
            print(f"✓ Variant import job started!")
            print(f"  Job ID: {job_id}")
            print(f"  Destination: {response['destinationName']}")
            print(f"  Run ID: {response.get('runId', 'N/A')}")
            
            # Monitor import
            print("\nMonitoring import progress...")
            self._monitor_variant_import(job_id)
            
            return job_id
            
        except ClientError as e:
            print(f"✗ Error importing VCF: {e}")
            raise
    
    def import_annotations(self,
                          s3_annotation_uri,
                          annotation_store_id,
                          role_arn,
                          format_to_header=None):
        """
        Import annotations into Annotation Store
        
        Args:
            s3_annotation_uri (str): S3 URI of annotation file (TSV or VCF)
            annotation_store_id (str): ID of annotation store
            role_arn (str): IAM role ARN for import job
            format_to_header (dict): Optional header mappings
            
        Returns:
            str: Import job ID
        """
        print("="*60)
        print("Importing Annotations")
        print("="*60)
        print(f"Source: {s3_annotation_uri}")
        print(f"Destination Store: {annotation_store_id}\n")
        
        try:
            # Prepare import request
            import_params = {
                'destinationName': annotation_store_id,
                'roleArn': role_arn,
                'items': [{
                    'source': s3_annotation_uri
                }]
            }
            
            if format_to_header:
                import_params['formatOptions'] = {
                    'tsvOptions': {
                        'readOptions': format_to_header
                    }
                }
            
            # Start import job
            response = self.omics_client.start_annotation_import_job(**import_params)
            
            job_id = response['id']
            print(f"✓ Annotation import job started!")
            print(f"  Job ID: {job_id}")
            print(f"  Destination: {response['destinationName']}")
            
            # Monitor import
            print("\nMonitoring import progress...")
            self._monitor_annotation_import(job_id)
            
            return job_id
            
        except ClientError as e:
            print(f"✗ Error importing annotations: {e}")
            raise
    
    def _monitor_variant_import(self, job_id, max_attempts=120):
        """Monitor variant import job progress"""
        for attempt in range(max_attempts):
            try:
                response = self.omics_client.get_variant_import_job(jobId=job_id)
                
                status = response['status']
                
                # Show progress
                if 'items' in response:
                    total = len(response['items'])
                    completed = sum(1 for item in response['items'] 
                                  if item.get('status') == 'COMPLETED')
                    progress = f" ({completed}/{total} items)"
                else:
                    progress = ""
                
                print(f"  [{attempt+1}/{max_attempts}] Status: {status}{progress}", end='\r')
                
                if status == 'COMPLETED':
                    print("\n\n✓ Variant import completed successfully!")
                    
                    # Show statistics
                    if 'runLeftNormalizationJobId' in response:
                        print(f"  Normalization completed")
                    
                    return
                    
                elif status == 'FAILED':
                    print("\n")
                    error_msg = response.get('statusMessage', 'Unknown error')
                    raise Exception(f"Import failed: {error_msg}")
                
                time.sleep(10)
                
            except ClientError as e:
                if attempt == max_attempts - 1:
                    raise
                time.sleep(10)
        
        raise Exception("Import job timed out")
    
    def _monitor_annotation_import(self, job_id, max_attempts=120):
        """Monitor annotation import job progress"""
        for attempt in range(max_attempts):
            try:
                response = self.omics_client.get_annotation_import_job(jobId=job_id)
                
                status = response['status']
                
                # Show progress
                if 'items' in response:
                    total = len(response['items'])
                    completed = sum(1 for item in response['items'] 
                                  if item.get('status') == 'COMPLETED')
                    progress = f" ({completed}/{total} items)"
                else:
                    progress = ""
                
                print(f"  [{attempt+1}/{max_attempts}] Status: {status}{progress}", end='\r')
                
                if status == 'COMPLETED':
                    print("\n\n✓ Annotation import completed successfully!")
                    return
                    
                elif status == 'FAILED':
                    print("\n")
                    error_msg = response.get('statusMessage', 'Unknown error')
                    raise Exception(f"Import failed: {error_msg}")
                
                time.sleep(10)
                
            except ClientError as e:
                if attempt == max_attempts - 1:
                    raise
                time.sleep(10)
        
        raise Exception("Import job timed out")
    
    def _wait_for_variant_store(self, store_name, max_attempts=60):
        """Wait for variant store to become active"""
        for attempt in range(max_attempts):
            try:
                response = self.omics_client.get_variant_store(name=store_name)
                status = response.get('status', 'UNKNOWN')
                
                if status == 'ACTIVE':
                    return
                elif status in ['FAILED', 'DELETING']:
                    raise Exception(f"Store entered {status} state")
                
                time.sleep(10)
            except ClientError as e:
                if attempt == max_attempts - 1:
                    raise
                time.sleep(10)
    
    def list_variant_stores(self):
        """List all variant stores"""
        print("\n" + "="*60)
        print("Variant Stores")
        print("="*60 + "\n")
        
        try:
            response = self.omics_client.list_variant_stores()
            
            stores = response.get('variantStores', [])
            if stores:
                for store in stores:
                    print(f"Store: {store['name']}")
                    print(f"  ID: {store['id']}")
                    print(f"  Status: {store.get('status', 'N/A')}")
                    print(f"  Reference: {store.get('reference', {}).get('referenceArn', 'N/A')}")
                    print(f"  Created: {store.get('creationTime', 'N/A')}")
                    print()
            else:
                print("No variant stores found")
                
        except ClientError as e:
            print(f"Error listing variant stores: {e}")
    
    def list_annotation_stores(self):
        """List all annotation stores"""
        print("\n" + "="*60)
        print("Annotation Stores")
        print("="*60 + "\n")
        
        try:
            response = self.omics_client.list_annotation_stores()
            
            stores = response.get('annotationStores', [])
            if stores:
                for store in stores:
                    print(f"Store: {store['name']}")
                    print(f"  ID: {store['id']}")
                    print(f"  Status: {store.get('status', 'N/A')}")
                    print(f"  Format: {store.get('storeFormat', 'N/A')}")
                    print(f"  Created: {store.get('creationTime', 'N/A')}")
                    print()
            else:
                print("No annotation stores found")
                
        except ClientError as e:
            print(f"Error listing annotation stores: {e}")
    
    def get_variant_store_stats(self, variant_store_name):
        """
        Get statistics about a variant store
        
        Args:
            variant_store_name (str): Name of the variant store
        """
        print("="*60)
        print(f"Variant Store Statistics: {variant_store_name}")
        print("="*60 + "\n")
        
        try:
            response = self.omics_client.get_variant_store(name=variant_store_name)
            
            print(f"Status: {response.get('status', 'N/A')}")
            print(f"Storage Size: {response.get('statusMessage', 'N/A')}")
            print(f"Reference: {response.get('reference', {}).get('referenceArn', 'N/A')}")
            
            # Additional stats if available
            if 'sseConfig' in response:
                print(f"Encryption: {response['sseConfig'].get('type', 'N/A')}")
            
        except ClientError as e:
            print(f"Error getting variant store stats: {e}")
    
    def _save_config(self):
        """Save updated configuration"""
        with open('/home/claude/healthomics_config.json', 'w') as f:
            json.dump(self.config, f, indent=2)


def create_sample_annotation_tsv():
    """
    Create a sample annotation TSV file for demonstration
    """
    print("Creating sample annotation TSV file...")
    
    annotation_data = """chr\tpos\tref\talt\tgene\tconsequence\tclinical_significance\tpopulation_af
chr1\t12345\tA\tG\tBRCA1\tmissense_variant\tPathogenic\t0.0001
chr2\t67890\tC\tT\tTP53\tstop_gained\tPathogenic\t0.0002
chr3\t11111\tG\tA\tEGFR\tsynonymous_variant\tBenign\t0.15
chr4\t22222\tT\tC\tKRAS\tmissense_variant\tLikely_pathogenic\t0.005
"""
    
    with open('/home/claude/sample_annotations.tsv', 'w') as f:
        f.write(annotation_data)
    
    print("✓ Sample annotation file created: /home/claude/sample_annotations.tsv")


def demonstrate_variant_management():
    """
    Demonstrate variant and annotation management
    """
    print("="*60)
    print("AWS HealthOmics Workshop - Module 4")
    print("Variant and Annotation Management")
    print("="*60)
    print("\nThis module demonstrates:")
    print("  1. Creating Variant Stores with references")
    print("  2. Importing VCF files")
    print("  3. Importing annotations")
    print("  4. Managing variant data")
    print("\n" + "="*60 + "\n")
    
    # Initialize manager
    manager = HealthOmicsVariantManager()
    
    # Show current stores
    manager.list_variant_stores()
    manager.list_annotation_stores()
    
    # Create sample annotation file
    print("\n")
    create_sample_annotation_tsv()
    
    print("\n" + "="*60)
    print("Variant Management Overview Complete!")
    print("="*60)
    print("\nKey Concepts:")
    print("  • VCF Import: Load variant calls from workflows")
    print("  • Annotations: Enrich variants with external data")
    print("  • Normalization: Standardize variant representations")
    print("  • Querying: Use Athena for powerful queries")
    print("\nNext Steps:")
    print("  1. Import your VCF files from workflow outputs")
    print("  2. Load annotation databases (ClinVar, dbSNP)")
    print("  3. Query variants using Athena (Module 5)")
    print("\n" + "="*60)


if __name__ == '__main__':
    demonstrate_variant_management()