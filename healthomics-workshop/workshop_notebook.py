#!/usr/bin/env python3
"""
AWS HealthOmics Workshop - VS Code Version

"""

# ============================================================================
# CELL 1: Import Required Libraries
# ============================================================================
print("=" * 70)
print("AWS HealthOmics Workshop - Setup")
print("=" * 70)
print("\nImporting required libraries...")

from datetime import datetime
import time
import boto3
import botocore.exceptions

# Create AWS service clients
omics = boto3.client('omics')
iam = boto3.client('iam')
s3 = boto3.client('s3')
sts = boto3.client('sts')

print("✓ Libraries imported successfully")
print("✓ AWS clients created")

# ============================================================================
# CELL 2: Configure Global Variables
# ============================================================================
print("\n" + "=" * 70)
print("Configuring Workshop Environment")
print("=" * 70)

# Get AWS account information
AWS_ACCOUNT_ID = sts.get_caller_identity()['Account']
AWS_REGION = omics.meta.region_name

print(f"\n✓ AWS Account ID: {AWS_ACCOUNT_ID}")
print(f"✓ AWS Region: {AWS_REGION}")

# Check/Create IAM Role
try:
    OMICS_JOB_ROLE_ARN = iam.get_role(RoleName='OmicsUnifiedJobRole')['Role']['Arn']
    print(f"✓ Found existing IAM role: OmicsUnifiedJobRole")
except iam.exceptions.NoSuchEntityException:
    print("⚠ IAM role 'OmicsUnifiedJobRole' not found")
    print("  This role is needed for HealthOmics operations")
    print("\n  Creating a basic IAM role for the workshop...")
    
    # Create the role (simplified version)
    try:
        trust_policy = {
            "Version": "2012-10-17",
            "Statement": [
                {
                    "Effect": "Allow",
                    "Principal": {
                        "Service": "omics.amazonaws.com"
                    },
                    "Action": "sts:AssumeRole"
                }
            ]
        }
        
        role_response = iam.create_role(
            RoleName='OmicsUnifiedJobRole',
            AssumeRolePolicyDocument=str(trust_policy),
            Description='Role for AWS HealthOmics workshop operations'
        )
        
        # Attach necessary policies
        iam.attach_role_policy(
            RoleName='OmicsUnifiedJobRole',
            PolicyArn='arn:aws:iam::aws:policy/AmazonS3FullAccess'
        )
        
        OMICS_JOB_ROLE_ARN = role_response['Role']['Arn']
        print(f"✓ Created IAM role: {OMICS_JOB_ROLE_ARN}")
        print("  Waiting 10 seconds for role to propagate...")
        time.sleep(10)
        
    except Exception as e:
        print(f"✗ Error creating IAM role: {e}")
        print("\n  Please create the role manually or ask your AWS administrator")
        print("  Required permissions: S3 access, HealthOmics access")
        OMICS_JOB_ROLE_ARN = f"arn:aws:iam::{AWS_ACCOUNT_ID}:role/OmicsUnifiedJobRole"

# Set bucket names
OMICS_WORKSHOP_BUCKET = f'aws-genomics-static-{AWS_REGION}/omics-workshop'
OMICS_OUTPUT_BUCKET = f'omics-output-{AWS_REGION}-{AWS_ACCOUNT_ID}'

print(f"\n✓ Workshop data bucket: s3://{OMICS_WORKSHOP_BUCKET}")
print(f"✓ Output bucket: s3://{OMICS_OUTPUT_BUCKET}")

# Create output bucket if it doesn't exist
print("\nChecking output bucket...")
try:
    s3.head_bucket(Bucket=OMICS_OUTPUT_BUCKET.split('/')[0])
    print(f"✓ Output bucket exists")
except:
    try:
        print(f"  Creating output bucket: {OMICS_OUTPUT_BUCKET.split('/')[0]}")
        if AWS_REGION == 'us-east-1':
            s3.create_bucket(Bucket=OMICS_OUTPUT_BUCKET.split('/')[0])
        else:
            s3.create_bucket(
                Bucket=OMICS_OUTPUT_BUCKET.split('/')[0],
                CreateBucketConfiguration={'LocationConstraint': AWS_REGION}
            )
        print(f"✓ Output bucket created")
    except Exception as e:
        print(f"⚠ Could not create bucket: {e}")
        print("  You may need to create it manually or use a different name")

# Generate timestamp for unique naming
TIMESTAMP = datetime.now().strftime('%Y%m%dT%H%M%S')
print(f"\n✓ Timestamp for this session: {TIMESTAMP}")

# ============================================================================
# Summary of Configuration
# ============================================================================
print("\n" + "=" * 70)
print("Configuration Summary")
print("=" * 70)
print(f"""
AWS Account ID:     {AWS_ACCOUNT_ID}
AWS Region:         {AWS_REGION}
IAM Role ARN:       {OMICS_JOB_ROLE_ARN}
Workshop Bucket:    s3://{OMICS_WORKSHOP_BUCKET}
Output Bucket:      s3://{OMICS_OUTPUT_BUCKET}
Session Timestamp:  {TIMESTAMP}
""")

# Save configuration for other scripts
import json

config = {
    'aws_account_id': AWS_ACCOUNT_ID,
    'aws_region': AWS_REGION,
    'omics_job_role_arn': OMICS_JOB_ROLE_ARN,
    'workshop_bucket': OMICS_WORKSHOP_BUCKET,
    'output_bucket': OMICS_OUTPUT_BUCKET,
    'timestamp': TIMESTAMP
}

with open('workshop_config.json', 'w') as f:
    json.dump(config, f, indent=2)

print("✓ Configuration saved to: workshop_config.json")

print("\n" + "=" * 70)
print("Setup Complete! Ready to proceed with workshop modules")
print("=" * 70)

# ============================================================================
# READY TO CONTINUE
# ============================================================================
print("""
Next Steps:
1. The environment is now configured
2. You can continue with the workshop modules
3. All variables are set and ready to use

Workshop Modules Available:
- Module 1: Create HealthOmics Stores (01_create_stores.py)
- Module 2: Upload Data (02_upload_data.py)
- Module 3: Run Workflows (03_workflow_management.py)
- Module 4: Manage Variants (04_variant_management.py)
- Module 5: Query with Athena (05_athena_queries.py)
""")