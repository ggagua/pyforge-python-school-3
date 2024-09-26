# Deploy to EC2 Instance

This repository contains a GitHub Actions workflow designed to automate the deployment of your application to an Amazon EC2 instance whenever changes are pushed to specific branches. This setup helps streamline the deployment process, ensuring that your application is always up-to-date with the latest code changes.

## Workflow Overview

The deployment workflow is triggered by:

- **Push Events**: When changes are pushed.
- **Pull Requests**: When pull requests are made.

### What the Workflow Does

1. **Checkout Repository**: 
   - The workflow starts by checking out the latest version of your repository, ensuring that the deployment has access to the most recent code.

2. **SSH into EC2 Instance**:
   - The workflow uses the `appleboy/ssh-action` to connect to your EC2 instance via SSH. It requires the instance's public IP, SSH username, and private SSH key, which should be securely stored as GitHub Secrets.
   - Upon successful connection, the workflow:
     - Updates the package manager on the EC2 instance to ensure all installed packages are current.
     - Installs Git if it's not already installed.
     - Checks for the existence of a project directory (`/home/<username>/myproject`):
       - If the directory **does not exist**, it creates the directory and clones the repository from GitHub.
       - If the directory **exists**, it navigates to the directory and pulls the latest changes from the `main` branch.
   - The workflow logs the status of the repository, confirming whether it is up-to-date with the latest changes.

### Workflow File

The workflow is defined in the `.github/workflows/deploy.yml` file and includes the following code:

```yaml
name: Deploy to EC2 Instance

on:
  push:
    branches: [main, ci]
  pull_request:
    branches: [main, ci]

jobs:
  deploy:
    runs-on: ubuntu-latest

    steps:
    - name: Checkout Repository
      uses: actions/checkout@v4

    - name: SSH into EC2 Instance
      uses: appleboy/ssh-action@v0.1.5
      with:
        host: ${{ secrets.EC2_HOST }}
        username: ${{ secrets.EC2_USER }}
        key: ${{ secrets.EC2_KEY }}
        port: 22
        debug: true
        script: |
          echo "Deploying to EC2 Instance..."
          
          sudo apt-get update
          sudo apt-get install -y git
          
          PROJECT_DIR="/home/${{ secrets.EC2_USER }}/myproject"

          if [ ! -d "$PROJECT_DIR" ]; then
            echo "Directory not found. Creating it and cloning repository..."
            mkdir -p "$PROJECT_DIR"
            cd "$PROJECT_DIR"
            git clone https://github.com/ggagua/pyforge-python-school-3.git . || echo "Failed to clone repository"
          else
            cd "$PROJECT_DIR"
            git pull origin main || echo "Failed to pull latest changes"
          fi

          git status
          echo "Repository is up to date."