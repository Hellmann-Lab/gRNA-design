# Read the requirements.txt file
requirements <- readLines("requirements.txt")

# Install each package
for (package in requirements) {
  install.packages(package)
}