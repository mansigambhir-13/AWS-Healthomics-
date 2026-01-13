#!/usr/bin/env python3
"""
Module 1: Create AWS HealthOmics Stores
This script demonstrates how to create all four types of HealthOmics stores using boto3
"""

import boto3
import time
import json
from botocore.exceptions import ClientError

class HealthOmicsStoreManager:
    """
    Manages creation and configuration of AWS HealthOmics stores
    
    This class provides methods to create and manage:
    - Sequence Stores: For raw sequencing data (FASTQ, BAM, CRAM)
    - Reference Stores: For reference genomes
    - Variant Stores: For variant calls (VCF files)
    - Annotation Stores: For variant annotations
    """
    
    def __init__(self, region='us-east-1', workshop_name='healthomics-workshop'):
        """
        Initialize the HealthOmics client
        
        Args:
            region (str): AWS region where stores will be created
            workshop_name (str): Prefix for store names
        """
        self.region = region
        self.workshop_name = workshop_name
        self.omics_client = boto3.client('omics', region_name=region)
        self.sts_client = boto3.client('sts', region_name=region)
        self.account_id = self.sts_client.get_caller_identity()['Account']
        
        print(f"Initialized HealthOmics client in region: {region}")
        print(f"AWS Account ID: {self.account_id}")
        print(f"Workshop name: {workshop_name}\n")
    
    def create_sequence_store(self):
        """
        Create a Sequence Store for raw sequencing data
        
        Sequence Stores automatically:
        - Compress data (50-70% reduction)
        - Encrypt data at rest
        - Enable fast retrieval
        
        Returns:
            str: Sequence Store ID
        """
        print("="*60)
        print("Creating Sequence Store")
        print("="*60)
        
        try:
            # Create the sequence store
            response = self.omics_client.create_sequence_store(
                name=f"{self.workshop_name}-sequence-store",
                description="Workshop sequence store for raw sequencing data",
                tags={
                    'Project': self.workshop_name,
                    'Purpose': 'Workshop',
                    'CreatedBy': 'HealthOmics-Workshop-Script'
                }
            )
            
            store_id = response['id']
            print(f"✓ Sequence Store created successfully!")
            print(f"  Store ID: {store_id}")
            print(f"  ARN: {response['arn']}")
            print(f"  Creation Time: {response['creationTime']}")
            
            # Wait for store to become active
            print("\nWaiting for Sequence Store to become active...")
            self._wait_for_sequence_store(store_id)
            print("✓ Sequence Store is now ACTIVE and ready to use")
            
            # Get detailed information
            store_info = self.get_sequence_store_info(store_id)
            self._print_store_details("Sequence Store", store_info)
            
            return store_id
            
        except ClientError as e:
            print(f"✗ Error creating Sequence Store: {e}")
            raise
    
    def create_reference_store(self):
        """
        Create a Reference Store for reference genomes
        
        Reference Stores:
        - Store reference sequences (e.g., GRCh38, GRCh37)
        - Support multiple reference versions
        - Enable efficient storage and retrieval
        
        Returns:
            str: Reference Store ID
        """
        print("\n" + "="*60)
        print("Creating Reference Store")
        print("="*60)
        
        try:
            response = self.omics_client.create_reference_store(
                name=f"{self.workshop_name}-reference-store",
                description="Workshop reference store for reference genomes",
                tags={
                    'Project': self.workshop_name,
                    'Purpose': 'Workshop',
                    'CreatedBy': 'HealthOmics-Workshop-Script'
                }
            )
            
            store_id = response['id']
            print(f"✓ Reference Store created successfully!")
            print(f"  Store ID: {store_id}")
            print(f"  ARN: {response['arn']}")
            print(f"  Creation Time: {response['creationTime']}")
            
            # Wait for store to become active
            print("\nWaiting for Reference Store to become active...")
            self._wait_for_reference_store(store_id)
            print("✓ Reference Store is now ACTIVE and ready to use")
            
            # Get detailed information
            store_info = self.get_reference_store_info(store_id)
            self._print_store_details("Reference Store", store_info)
            
            return store_id
            
        except ClientError as e:
            print(f"✗ Error creating Reference Store: {e}")
            raise
    
    def create_annotation_store(self, store_format='TSV'):
        """
        Create an Annotation Store for variant annotations
        
        Annotation Stores:
        - Store variant annotations (e.g., ClinVar, dbSNP)
        - Support TSV and VCF formats
        - Enable complex queries
        
        Args:
            store_format (str): Format of annotations (TSV or VCF)
            
        Returns:
            str: Annotation Store ID
        """
        print("\n" + "="*60)
        print("Creating Annotation Store")
        print("="*60)
        
        try:
            # Define schema for TSV format
            schema = [
                {
                    'name': 'chr',
                    'type': 'STRING',
                    'description': 'Chromosome'
                },
                {
                    'name': 'pos',
                    'type': 'LONG',
                    'description': 'Position'
                },
                {
                    'name': 'ref',
                    'type': 'STRING',
                    'description': 'Reference allele'
                },
                {
                    'name': 'alt',
                    'type': 'STRING',
                    'description': 'Alternate allele'
                },
                {
                    'name': 'gene',
                    'type': 'STRING',
                    'description': 'Gene name'
                },
                {
                    'name': 'consequence',
                    'type': 'STRING',
                    'description': 'Variant consequence'
                }
            ]
            
            response = self.omics_client.create_annotation_store(
                name=f"{self.workshop_name}-annotation-store",
                description="Workshop annotation store for variant annotations",
                storeFormat=store_format,
                tags={
                    'Project': self.workshop_name,
                    'Purpose': 'Workshop',
                    'CreatedBy': 'HealthOmics-Workshop-Script'
                }
            )
            
            store_id = response['id']
            print(f"✓ Annotation Store created successfully!")
            print(f"  Store ID: {store_id}")
            print(f"  Store Name: {response['name']}")
            print(f"  Format: {store_format}")
            print(f"  Creation Time: {response['creationTime']}")
            
            # Wait for store to become active
            print("\nWaiting for Annotation Store to become active...")
            self._wait_for_annotation_store(response['name'])
            print("✓ Annotation Store is now ACTIVE and ready to use")
            
            return store_id
            
        except ClientError as e:
            print(f"✗ Error creating Annotation Store: {e}")
            raise
    
    def _wait_for_sequence_store(self, store_id, max_attempts=30):
        """Wait for sequence store to become active"""
        for attempt in range(max_attempts):
            try:
                response = self.omics_client.get_sequence_store(id=store_id)
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
    
    def _wait_for_reference_store(self, store_id, max_attempts=30):
        """Wait for reference store to become active"""
        for attempt in range(max_attempts):
            try:
                response = self.omics_client.get_reference_store(id=store_id)
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
    
    def _wait_for_annotation_store(self, store_name, max_attempts=30):
        """Wait for annotation store to become active"""
        for attempt in range(max_attempts):
            try:
                response = self.omics_client.get_annotation_store(name=store_name)
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
    
    def get_sequence_store_info(self, store_id):
        """Get detailed information about a sequence store"""
        return self.omics_client.get_sequence_store(id=store_id)
    
    def get_reference_store_info(self, store_id):
        """Get detailed information about a reference store"""
        return self.omics_client.get_reference_store(id=store_id)
    
    def _print_store_details(self, store_type, info):
        """Print formatted store details"""
        print(f"\n{store_type} Details:")
        print(f"  Status: {info.get('status', 'N/A')}")
        print(f"  ARN: {info.get('arn', 'N/A')}")
        if 'sseConfig' in info:
            print(f"  Encryption: {info['sseConfig'].get('type', 'N/A')}")
    
    def list_all_stores(self):
        """List all stores created in this workshop"""
        print("\n" + "="*60)
        print("Summary of All Stores")
        print("="*60)
        
        # List Sequence Stores
        print("\nSequence Stores:")
        try:
            response = self.omics_client.list_sequence_stores()
            stores = [s for s in response.get('sequenceStores', []) 
                     if self.workshop_name in s.get('name', '')]
            if stores:
                for store in stores:
                    print(f"  • {store['name']}")
                    print(f"    ID: {store['id']}")
                    print(f"    Status: {store.get('status', 'N/A')}")
            else:
                print("  No sequence stores found")
        except ClientError as e:
            print(f"  Error listing: {e}")
        
        # List Reference Stores
        print("\nReference Stores:")
        try:
            response = self.omics_client.list_reference_stores()
            stores = [s for s in response.get('referenceStores', []) 
                     if self.workshop_name in s.get('name', '')]
            if stores:
                for store in stores:
                    print(f"  • {store['name']}")
                    print(f"    ID: {store['id']}")
                    print(f"    Status: {store.get('status', 'N/A')}")
            else:
                print("  No reference stores found")
        except ClientError as e:
            print(f"  Error listing: {e}")
        
        # List Annotation Stores
        print("\nAnnotation Stores:")
        try:
            response = self.omics_client.list_annotation_stores()
            stores = [s for s in response.get('annotationStores', []) 
                     if self.workshop_name in s.get('name', '')]
            if stores:
                for store in stores:
                    print(f"  • {store['name']}")
                    print(f"    ID: {store['id']}")
                    print(f"    Status: {store.get('status', 'N/A')}")
                    print(f"    Format: {store.get('storeFormat', 'N/A')}")
            else:
                print("  No annotation stores found")
        except ClientError as e:
            print(f"  Error listing: {e}")
    
    def save_store_ids(self, sequence_id=None, reference_id=None, annotation_id=None):
        """Save store IDs to a configuration file"""
        config = {
            'region': self.region,
            'account_id': self.account_id,
            'workshop_name': self.workshop_name,
            'sequence_store_id': sequence_id,
            'reference_store_id': reference_id,
            'annotation_store_id': annotation_id
        }
        
        with open('/home/claude/healthomics_config.json', 'w') as f:
            json.dump(config, f, indent=2)
        
        print("\n✓ Store IDs saved to: /home/claude/healthomics_config.json")


def main():
    """
    Main function to create all HealthOmics stores
    """
    print("="*60)
    print("AWS HealthOmics Workshop - Module 1")
    print("Creating HealthOmics Stores")
    print("="*60)
    print("\nThis script will create:")
    print("  1. Sequence Store - for raw sequencing data")
    print("  2. Reference Store - for reference genomes")
    print("  3. Annotation Store - for variant annotations")
    print("\n" + "="*60 + "\n")
    
    # Initialize manager
    manager = HealthOmicsStoreManager(region='us-east-1')
    
    try:
        # Create stores
        sequence_id = manager.create_sequence_store()
        reference_id = manager.create_reference_store()
        annotation_id = manager.create_annotation_store()
        
        # Save configuration
        manager.save_store_ids(
            sequence_id=sequence_id,
            reference_id=reference_id,
            annotation_id=annotation_id
        )
        
        # Display summary
        manager.list_all_stores()
        
        print("\n" + "="*60)
        print("✓ Store Creation Complete!")
        print("="*60)
        print("\nNext Steps:")
        print("  1. Upload a reference genome (Module 2a)")
        print("  2. Upload sequencing data (Module 2b)")
        print("  3. Run a workflow (Module 3)")
        print("  4. Load and query variants (Module 4-5)")
        print("\n" + "="*60)
        
    except Exception as e:
        print(f"\n✗ Error during store creation: {e}")
        print("Please check your AWS credentials and permissions.")
        raise


if __name__ == '__main__':
    main()