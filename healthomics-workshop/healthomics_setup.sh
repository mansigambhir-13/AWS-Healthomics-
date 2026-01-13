#!/bin/bash

# AWS HealthOmics Workshop - Environment Setup Script
# This script sets up your environment for the workshop

echo "=========================================="
echo "AWS HealthOmics Workshop Setup"
echo "=========================================="

# Set variables
export AWS_REGION="us-east-1"  # Change to your region if different
export WORKSHOP_NAME="healthomics-workshop"

# Create a directory structure for the workshop
echo "Creating workshop directory structure..."
mkdir -p ~/healthomics-workshop/{data,scripts,workflows,outputs}

# Check AWS CLI installation
echo "Checking AWS CLI..."
if ! command -v aws &> /dev/null; then
    echo "AWS CLI not found. Please install it first."
    echo "Visit: https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html"
    exit 1
fi

# Check AWS credentials
echo "Checking AWS credentials..."
if ! aws sts get-caller-identity &> /dev/null; then
    echo "AWS credentials not configured. Please configure them first."
    echo "Run: aws configure"
    exit 1
fi

# Get AWS account ID
export AWS_ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)
echo "AWS Account ID: $AWS_ACCOUNT_ID"

# Check if HealthOmics is available in the region
echo "Checking HealthOmics availability in region: $AWS_REGION"
if aws omics list-sequence-stores --region $AWS_REGION &> /dev/null; then
    echo "✓ HealthOmics is available in $AWS_REGION"
else
    echo "✗ HealthOmics may not be available in $AWS_REGION"
    echo "Available regions: us-east-1, us-west-2, eu-west-1, eu-west-2, eu-central-1, ap-southeast-1"
fi

# Export environment variables to a file
cat > ~/healthomics-workshop/env_vars.sh << EOF
export AWS_REGION="$AWS_REGION"
export AWS_ACCOUNT_ID="$AWS_ACCOUNT_ID"
export WORKSHOP_NAME="$WORKSHOP_NAME"
export WORKSHOP_DIR="$HOME/healthomics-workshop"
EOF

echo ""
echo "=========================================="
echo "Setup Complete!"
echo "=========================================="
echo "Workshop directory: ~/healthomics-workshop"
echo "To load environment variables, run:"
echo "source ~/healthomics-workshop/env_vars.sh"
echo ""
echo "Next steps:"
echo "1. Create HealthOmics stores"
echo "2. Upload sample data"
echo "3. Run workflows"
echo "=========================================="