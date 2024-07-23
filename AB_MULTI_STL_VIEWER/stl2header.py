import os

def file_to_cpp_header(input_filename, output_filename, array_name):
    try:
        # Read the binary content of the input file
        with open(input_filename, 'rb') as file:
            content = file.read()

        # Convert the binary content to a C++ unsigned char array
        array_content = ', '.join(f'0x{byte:02x}' for byte in content)

        # Prepare the header content
        header_content = f"""#pragma once

const unsigned char {array_name}[] = {{
    {array_content}
}};

"""

        # Write the header file
        with open(output_filename, 'w') as file:
            file.write(header_content)
        
        print(f"Successfully converted {input_filename} to {output_filename}")

    except IOError as e:
        print(f"Error processing file: {e}")

def convert_stl_to_headers(directory):
    # Iterate through all files in the given directory
    for filename in os.listdir(directory):
        if filename.lower().endswith('.stl'):
            input_path = os.path.join(directory, filename)
            output_filename = filename.replace('.stl', '_stl.h')
            output_path = os.path.join(directory, output_filename)
            array_name = os.path.splitext(filename)[0]
            file_to_cpp_header(input_path, output_path, array_name)

# Get the current script directory
current_directory = os.path.dirname(os.path.abspath(__file__))

# Convert all STL files in the script's directory
convert_stl_to_headers(current_directory)
