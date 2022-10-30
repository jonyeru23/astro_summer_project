from get_model import *

def main():
    sheet_name = input('Enter sheet name:')
    steps = int(input('Enter steps:'))
    sampler = Sampler('combined_data.xlsx', sheet_name=sheet_name)
    sampler.write_sampler(f"{sheet_name}.h5", steps=steps)


if __name__ == '__main__':
    main()

