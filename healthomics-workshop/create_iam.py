#!/usr/bin/env python3
"""
Create IAM Role for AWS HealthOmics Workshop
This script creates the OmicsUnifiedJobRole needed for the workshop
"""

import boto3
import json
import time
from botocore.exceptions import ClientError

def create_omics_role():
    """
    Create the OmicsUnifiedJobRole with necessary permissions
    """
    iam = boto3.client('iam')
    
    print("="*70)
    print("Creating IAM Role for AWS HealthOmics Workshop")
    print("="*70)
    
    role_name = 'OmicsUnifiedJobRole'
    
    # Check if role already exists
    try:
        existing_role = iam.get_role(RoleName=role_name)
        print(f"\n✓ Role '{role_name}' already exists")
        print(f"  ARN: {existing_role['Role']['Arn']}")
        print("\nNo action needed. You can proceed with the workshop.")
        return existing_role['Role']['Arn']
    except iam.exceptions.NoSuchEntityException:
        print(f"\nRole '{role_name}' not found. Creating it now...")
    
    # Trust policy - allows HealthOmics service to assume this role
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
    
    try:
        # Create the role
        print("\nStep 1: Creating IAM role...")
        role_response = iam.create_role(
            RoleName=role_name,
            AssumeRolePolicyDocument=json.dumps(trust_policy),
            Description='IAM role for AWS HealthOmics workshop operations',
            Tags=[
                {'Key': 'Purpose', 'Value': 'HealthOmics-Workshop'},
                {'Key': 'ManagedBy', 'Value': 'Workshop-Script'}
            ]
        )
        
        role_arn = role_response['Role']['Arn']
        print(f"✓ Role created: {role_arn}")
        
        # Attach managed policies
        print("\nStep 2: Attaching policies...")
        
        # S3 access for reading/writing data
        print("  - Attaching S3 access policy...")
        iam.attach_role_policy(
            RoleName=role_name,
            PolicyArn='arn:aws:iam::aws:policy/AmazonS3FullAccess'
        )
        print("    ✓ S3 access granted")
        
        # Create custom policy for HealthOmics operations
        print("  - Creating custom HealthOmics policy...")
        
        omics_policy = {
            "Version": "2012-10-17",
            "Statement": [
                {
                    "Effect": "Allow",
                    "Action": [
                        "omics:*"
                    ],
                    "Resource": "*"
                },
                {
                    "Effect": "Allow",
                    "Action": [
                        "logs:CreateLogGroup",
                        "logs:CreateLogStream",
                        "logs:PutLogEvents"
                    ],
                    "Resource": "arn:aws:logs:*:*:*"
                },
                {
                    "Effect": "Allow",
                    "Action": [
                        "ecr:GetAuthorizationToken",
                        "ecr:BatchCheckLayerAvailability",
                        "ecr:GetDownloadUrlForLayer",
                        "ecr:BatchGetImage"
                    ],
                    "Resource": "*"
                }
            ]
        }
        
        policy_name = 'OmicsWorkshopPolicy'
        
        try:
            policy_response = iam.create_policy(
                PolicyName=policy_name,
                PolicyDocument=json.dumps(omics_policy),
                Description='Policy for HealthOmics workshop operations'
            )
            policy_arn = policy_response['Policy']['Arn']
            print(f"    ✓ Policy created: {policy_arn}")
            
            # Attach the custom policy
            iam.attach_role_policy(
                RoleName=role_name,
                PolicyArn=policy_arn
            )
            print("    ✓ Custom policy attached")
            
        except ClientError as e:
            if e.response['Error']['Code'] == 'EntityAlreadyExists':
                # Policy already exists, get its ARN
                account_id = boto3.client('sts').get_caller_identity()['Account']
                policy_arn = f"arn:aws:iam::{account_id}:policy/{policy_name}"
                print(f"    ✓ Policy already exists: {policy_arn}")
                
                # Attach it
                iam.attach_role_policy(
                    RoleName=role_name,
                    PolicyArn=policy_arn
                )
                print("    ✓ Existing policy attached")
            else:
                raise
        
        print("\nStep 3: Waiting for role to propagate...")
        print("  (This takes a few seconds...)")
        time.sleep(10)
        print("  ✓ Role should now be ready")
        
        print("\n" + "="*70)
        print("✓ IAM Role Setup Complete!")
        print("="*70)
        print(f"\nRole Name: {role_name}")
        print(f"Role ARN: {role_arn}")
        print("\nThis role can now be used for:")
        print("  • Reading/writing S3 data")
        print("  • HealthOmics operations")
        print("  • Workflow execution")
        print("  • Logging to CloudWatch")
        
        print("\nYou can now proceed with the workshop!")
        print("="*70)
        
        return role_arn
        
    except ClientError as e:
        error_code = e.response['Error']['Code']
        error_msg = e.response['Error']['Message']
        
        print(f"\n✗ Error creating IAM role: {error_code}")
        print(f"  {error_msg}")
        
        if error_code == 'AccessDenied':
            print("\n⚠ You don't have permission to create IAM roles.")
            print("  Please ask your AWS administrator to:")
            print("  1. Create a role named 'OmicsUnifiedJobRole'")
            print("  2. Attach AmazonS3FullAccess policy")
            print("  3. Add HealthOmics permissions")
            print("  4. Allow omics.amazonaws.com to assume the role")
        
        raise

def verify_role():
    """Verify the role was created successfully"""
    iam = boto3.client('iam')
    role_name = 'OmicsUnifiedJobRole'
    
    print("\n" + "="*70)
    print("Verifying Role Configuration")
    print("="*70)
    
    try:
        # Get role details
        role = iam.get_role(RoleName=role_name)
        print(f"\n✓ Role exists: {role_name}")
        print(f"  ARN: {role['Role']['Arn']}")
        print(f"  Created: {role['Role']['CreateDate']}")
        
        # List attached policies
        print("\nAttached Policies:")
        policies = iam.list_attached_role_policies(RoleName=role_name)
        for policy in policies['AttachedPolicies']:
            print(f"  • {policy['PolicyName']}")
        
        print("\n✓ Role is properly configured")
        print("="*70)
        
    except ClientError as e:
        print(f"\n✗ Error verifying role: {e}")


def main():
    """Main execution"""
    try:
        # Check AWS credentials
        sts = boto3.client('sts')
        identity = sts.get_caller_identity()
        print(f"\nAWS Account: {identity['Account']}")
        print(f"User/Role: {identity['Arn']}")
        
        # Create the role
        role_arn = create_omics_role()
        
        # Verify it
        verify_role()
        
        # Save to config
        print("\nSaving role ARN to configuration...")
        try:
            with open('workshop_config.json', 'r') as f:
                config = json.load(f)
            config['omics_job_role_arn'] = role_arn
            with open('workshop_config.json', 'w') as f:
                json.dump(config, f, indent=2)
            print("✓ Configuration updated")
        except FileNotFoundError:
            # Create new config
            config = {
                'omics_job_role_arn': role_arn,
                'aws_account_id': identity['Account']
            }
            with open('workshop_config.json', 'w') as f:
                json.dump(config, f, indent=2)
            print("✓ Configuration created")
        
        print("\n" + "="*70)
        print("Setup complete! You can now run the workshop.")
        print("="*70)
        
    except ClientError as e:
        print(f"\n✗ AWS Error: {e}")
        print("\nPlease check:")
        print("  1. AWS credentials are configured (aws configure)")
        print("  2. You have IAM permissions to create roles")
        print("  3. You're using the correct AWS region")
        return 1
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        return 1
    
    return 0


if __name__ == '__main__':
    exit(main())