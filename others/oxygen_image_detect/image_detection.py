import csv
from PIL import Image

# YMC2007 data ---------------------------------------------------------------------

def extract_blue_pixels(image_path):
    img = Image.open(image_path)
    pixels = img.load()
    width, height = img.size
    blue_pixels = []

    for y in range(height):
        for x in range(width):
            r, g, b = pixels[x, y]
            if 30 <= r <= 105 and 55 <= g <= 120 and 70 <= b <= 140 and b > r and b > g:
                blue_pixels.append((x, y, r, g, b))

    return blue_pixels

image_path = 'LC-MS-1.jpg'
blue_pixels = extract_blue_pixels(image_path)

csv_file = 'LC-MS_all.csv'
with open(csv_file, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['x', 'y', 'r', 'g', 'b']) 
    writer.writerows(blue_pixels)



image_path = 'TOFMS-1.jpg'
blue_pixels = extract_blue_pixels(image_path)

csv_file = 'TOFMS_all.csv'
with open(csv_file, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['x', 'y', 'r', 'g', 'b']) 
    writer.writerows(blue_pixels)



# YMC2014 data ---------------------------------------------------------------------------------

def extract_blue_pixels2(image_path):
    img = Image.open(image_path)
    pixels = img.load()
    width, height = img.size
    blue_pixels = []

    for y in range(height):
        for x in range(width):
            r, g, b = pixels[x, y]
            if r <= 50 and g <= 50 and 200 <= b and b > r and b > g:
                blue_pixels.append((x, y, r, g, b))

    return blue_pixels


image_path = 'YMC2014_oxygen_chipseq.jpg'
blue_pixels = extract_blue_pixels2(image_path)

csv_file = 'YMC2014_oxygen_chipseq_all.csv'
with open(csv_file, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['x', 'y', 'r', 'g', 'b']) 
    writer.writerows(blue_pixels)



image_path = 'YMC2014_oxygen_rnaseq.jpg'
blue_pixels = extract_blue_pixels2(image_path)

csv_file = 'YMC2014_oxygen_rnaseq_all.csv'
with open(csv_file, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['x', 'y', 'r', 'g', 'b'])  
    writer.writerows(blue_pixels)



# YMC2005 data ---------------------------------------------------------------------

def extract_blue_pixels3(image_path):
    img = Image.open(image_path)
    pixels = img.load()
    width, height = img.size
    blue_pixels = []

    for y in range(height):
        for x in range(width):
            r, g, b = pixels[x, y]
            if r >= 150 and g >= 200 and 200 <= b and b > r and g > r:
                blue_pixels.append((x, y, r, g, b))

    return blue_pixels

image_path = 'YMC2005.jpg'
blue_pixels = extract_blue_pixels3(image_path)

csv_file = 'YMC2005_oxygen_all.csv'
with open(csv_file, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['x', 'y', 'r', 'g', 'b'])
    writer.writerows(blue_pixels)



