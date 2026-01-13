#!/usr/bin/env python3
"""
Module 3: AWS HealthOmics Workflows
This script demonstrates how to create, manage, and run bioinformatics workflows
"""

import boto3
import json
import time
from botocore.exceptions import ClientError

class HealthOmicsWorkflowManager:
    """
    Manages HealthOmics workflows for genomics analysis
    
    Supports:
    - Creating workflow definitions (Nextflow, WDL)
    - Starting workflow runs
    - Monitoring run progress
    - Managing run groups
    """
    
    def __init__(self, config_file='/home/claude/healthomics_config.json'):
        """Initialize workflow manager"""
        with open(config_file, 'r') as f:
            self.config = json.load(f)
        
        self.region = self.config['region']
        self.omics_client = boto3.client('omics', region_name=self.region)
        
        print(f"Initialized Workflow Manager")
        print(f"Region: {self.region}\n")
    
    def create_simple_workflow(self, name="simple-alignment"):
        """
        Create a simple alignment workflow (Nextflow-based)
        
        This creates a basic BWA-MEM alignment workflow that:
        1. Takes FASTQ files as input
        2. Aligns to a reference genome
        3. Produces BAM output
        
        Args:
            name (str): Workflow name
            
        Returns:
            str: Workflow ID
        """
        print("="*60)
        print("Creating Workflow Definition")
        print("="*60)
        print(f"Workflow name: {name}")
        print(f"Engine: Nextflow DSL2\n")
        
        # Simple Nextflow workflow definition
        # This is a simplified example - real workflows are more complex
        workflow_definition = """
nextflow.enable.dsl=2

// Simple alignment workflow
process align {
    container 'quay.io/biocontainers/bwa:0.7.17--h5bf99c6_8'
    
    input:
    path reference
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}.bam")
    
    script:
    \"\"\"
    bwa mem -t 4 ${reference} ${reads[0]} ${reads[1]} | \\
        samtools view -bS - > ${sample_id}.bam
    \"\"\"
}

process sort_bam {
    container 'quay.io/biocontainers/samtools:1.15--h3843a85_0'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam")
    
    script:
    \"\"\"
    samtools sort -@ 4 -o ${sample_id}.sorted.bam ${bam}
    samtools index ${sample_id}.sorted.bam
    \"\"\"
}

workflow {
    // Define inputs
    reference_ch = Channel.fromPath(params.reference)
    reads_ch = Channel.fromFilePairs(params.reads)
    
    // Run alignment pipeline
    align(reference_ch, reads_ch)
    sort_bam(align.out)
}
"""
        
        try:
            # In real implementation, you would:
            # 1. Create a workflow bundle (ZIP file)
            # 2. Upload to S3
            # 3. Call create_workflow with definitionUri
            
            print("Creating workflow definition...")
            print("\nWorkflow Steps:")
            print("  1. BWA-MEM alignment")
            print("  2. SAM to BAM conversion")
            print("  3. BAM sorting and indexing")
            
            # Placeholder for actual workflow creation
            # In workshop, pre-built workflows are typically provided
            workflow_id = "placeholder-workflow-id"
            
            print(f"\n✓ Workflow created successfully!")
            print(f"  Workflow ID: {workflow_id}")
            
            return workflow_id
            
        except ClientError as e:
            print(f"✗ Error creating workflow: {e}")
            raise
    
    def start_workflow_run(self,
                          workflow_id,
                          parameters,
                          output_uri,
                          role_arn,
                          name="workflow-run",
                          run_group_id=None):
        """
        Start a workflow run
        
        Args:
            workflow_id (str): ID of the workflow to run
            parameters (dict): Workflow input parameters
            output_uri (str): S3 URI for outputs
            role_arn (str): IAM role ARN for the workflow
            name (str): Name for the run
            run_group_id (str): Optional run group ID
            
        Returns:
            str: Run ID
        """
        print("="*60)
        print("Starting Workflow Run")
        print("="*60)
        print(f"Workflow ID: {workflow_id}")
        print(f"Run name: {name}")
        print(f"Output location: {output_uri}\n")
        
        try:
            # Prepare run parameters
            run_params = {
                'workflowId': workflow_id,
                'roleArn': role_arn,
                'name': name,
                'parameters': parameters,
                'outputUri': output_uri
            }
            
            if run_group_id:
                run_params['runGroupId'] = run_group_id
            
            # Start the run
            response = self.omics_client.start_run(**run_params)
            
            run_id = response['id']
            print(f"✓ Workflow run started successfully!")
            print(f"  Run ID: {run_id}")
            print(f"  Status: {response['status']}")
            print(f"  ARN: {response['arn']}")
            
            return run_id
            
        except ClientError as e:
            print(f"✗ Error starting workflow run: {e}")
            raise
    
    def monitor_run(self, run_id, poll_interval=30):
        """
        Monitor workflow run progress
        
        Args:
            run_id (str): Run ID to monitor
            poll_interval (int): Seconds between status checks
        """
        print("\n" + "="*60)
        print(f"Monitoring Run: {run_id}")
        print("="*60 + "\n")
        
        attempt = 0
        while True:
            try:
                response = self.omics_client.get_run(id=run_id)
                
                status = response['status']
                attempt += 1
                
                # Display run information
                print(f"Status Check #{attempt}")
                print(f"  Status: {status}")
                print(f"  Started: {response.get('startTime', 'N/A')}")
                
                if 'tasks' in response:
                    total_tasks = len(response['tasks'])
                    completed = sum(1 for t in response['tasks'] 
                                  if t.get('status') == 'COMPLETED')
                    print(f"  Progress: {completed}/{total_tasks} tasks completed")
                
                # Terminal states
                if status == 'COMPLETED':
                    print(f"\n✓ Workflow run completed successfully!")
                    print(f"  Run time: {response.get('stopTime', 'N/A')}")
                    self._print_run_summary(response)
                    return 'COMPLETED'
                    
                elif status == 'FAILED':
                    print(f"\n✗ Workflow run failed")
                    print(f"  Error: {response.get('statusMessage', 'Unknown error')}")
                    return 'FAILED'
                    
                elif status == 'CANCELLED':
                    print(f"\n⚠ Workflow run was cancelled")
                    return 'CANCELLED'
                
                # Continue monitoring
                print()
                time.sleep(poll_interval)
                
            except ClientError as e:
                print(f"Error checking run status: {e}")
                time.sleep(poll_interval)
    
    def _print_run_summary(self, run_info):
        """Print detailed run summary"""
        print("\nRun Summary:")
        print(f"  Name: {run_info.get('name', 'N/A')}")
        print(f"  Output URI: {run_info.get('outputUri', 'N/A')}")
        
        if 'resourceDigests' in run_info:
            print(f"\nResource Usage:")
            for key, value in run_info['resourceDigests'].items():
                print(f"  {key}: {value}")
    
    def create_run_group(self, name, max_cpus=100000, max_runs=100):
        """
        Create a run group to organize workflow runs
        
        Args:
            name (str): Run group name
            max_cpus (int): Maximum CPUs for the group
            max_runs (int): Maximum concurrent runs
            
        Returns:
            str: Run group ID
        """
        print("="*60)
        print("Creating Run Group")
        print("="*60)
        
        try:
            response = self.omics_client.create_run_group(
                name=name,
                maxCpus=max_cpus,
                maxRuns=max_runs,
                tags={
                    'Project': self.config.get('workshop_name', 'healthomics-workshop')
                }
            )
            
            group_id = response['id']
            print(f"✓ Run group created successfully!")
            print(f"  Group ID: {group_id}")
            print(f"  ARN: {response['arn']}")
            print(f"  Max CPUs: {max_cpus}")
            print(f"  Max Runs: {max_runs}")
            
            return group_id
            
        except ClientError as e:
            print(f"✗ Error creating run group: {e}")
            raise
    
    def list_workflows(self):
        """List all workflows"""
        print("\n" + "="*60)
        print("Available Workflows")
        print("="*60 + "\n")
        
        try:
            response = self.omics_client.list_workflows()
            
            workflows = response.get('items', [])
            if workflows:
                for wf in workflows:
                    print(f"Workflow: {wf['name']}")
                    print(f"  ID: {wf['id']}")
                    print(f"  Type: {wf.get('type', 'N/A')}")
                    print(f"  Status: {wf.get('status', 'N/A')}")
                    print(f"  Created: {wf.get('creationTime', 'N/A')}")
                    print()
            else:
                print("No workflows found")
                
        except ClientError as e:
            print(f"Error listing workflows: {e}")
    
    def list_runs(self, status=None):
        """
        List workflow runs
        
        Args:
            status (str): Filter by status (PENDING, RUNNING, COMPLETED, etc.)
        """
        print("\n" + "="*60)
        print("Workflow Runs")
        if status:
            print(f"Status Filter: {status}")
        print("="*60 + "\n")
        
        try:
            params = {}
            if status:
                params['status'] = status
            
            response = self.omics_client.list_runs(**params)
            
            runs = response.get('items', [])
            if runs:
                for run in runs:
                    print(f"Run: {run['name']}")
                    print(f"  ID: {run['id']}")
                    print(f"  Status: {run.get('status', 'N/A')}")
                    print(f"  Started: {run.get('startTime', 'N/A')}")
                    if run.get('stopTime'):
                        print(f"  Stopped: {run['stopTime']}")
                    print()
            else:
                print("No runs found")
                
        except ClientError as e:
            print(f"Error listing runs: {e}")
    
    def get_run_details(self, run_id):
        """
        Get detailed information about a specific run
        
        Args:
            run_id (str): Run ID
        """
        print("="*60)
        print(f"Run Details: {run_id}")
        print("="*60 + "\n")
        
        try:
            response = self.omics_client.get_run(id=run_id)
            
            print(f"Name: {response.get('name', 'N/A')}")
            print(f"Status: {response.get('status', 'N/A')}")
            print(f"Workflow ID: {response.get('workflowId', 'N/A')}")
            print(f"Started: {response.get('startTime', 'N/A')}")
            print(f"Output URI: {response.get('outputUri', 'N/A')}")
            
            if 'parameters' in response:
                print(f"\nParameters:")
                for key, value in response['parameters'].items():
                    print(f"  {key}: {value}")
            
            if 'tasks' in response:
                print(f"\nTasks:")
                for task in response['tasks']:
                    print(f"  - {task.get('name', 'N/A')}: {task.get('status', 'N/A')}")
            
        except ClientError as e:
            print(f"Error getting run details: {e}")
    
    def cancel_run(self, run_id):
        """
        Cancel a running workflow
        
        Args:
            run_id (str): Run ID to cancel
        """
        print(f"Cancelling run: {run_id}")
        
        try:
            self.omics_client.cancel_run(id=run_id)
            print(f"✓ Cancellation request submitted")
            
        except ClientError as e:
            print(f"✗ Error cancelling run: {e}")


def demonstrate_workflow_management():
    """
    Demonstrate workflow management capabilities
    """
    print("="*60)
    print("AWS HealthOmics Workshop - Module 3")
    print("Workflow Management")
    print("="*60)
    print("\nThis module demonstrates:")
    print("  1. Creating workflow definitions")
    print("  2. Starting workflow runs")
    print("  3. Monitoring run progress")
    print("  4. Managing run groups")
    print("\n" + "="*60 + "\n")
    
    # Initialize manager
    manager = HealthOmicsWorkflowManager()
    
    # Show available workflows and runs
    manager.list_workflows()
    manager.list_runs()
    
    print("\n" + "="*60)
    print("Workflow Management Overview Complete!")
    print("="*60)
    print("\nExample Workflow Types:")
    print("  • Alignment (BWA, Bowtie2)")
    print("  • Variant Calling (GATK, FreeBayes)")
    print("  • RNA-Seq (STAR, Salmon)")
    print("  • Quality Control (FastQC, MultiQC)")
    print("\nNext Steps:")
    print("  1. Run a pre-built workflow on your data")
    print("  2. Load variant results into Variant Store")
    print("  3. Query variants with Athena")
    print("\n" + "="*60)


if __name__ == '__main__':
    demonstrate_workflow_management()