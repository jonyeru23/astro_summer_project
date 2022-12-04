from make_figs import *
from get_model import *
import os

"""
i want to get the *.h5 file and the graphs in the same run in the nove, it takes too long

the command is in this way because this way you can open several screen sessions in the nova and let it run all together
"""

def main():
    """
    this is a heavy model, it will propably take you a few days to run.
    remember to open different screen sessions to make it run parallel

    @pre: models in os.listdir()
    """
    sheet_name = input('Enter sheet name:')
    steps = int(input('Enter steps:'))

    data_path = "combined_data.xlsx"
    folder_path = os.path.join("models", f"{sheet_name} {steps}")
    model_path = os.path.join(folder_path, f"{sheet_name}_{steps}.h5")
    os.mkdir(folder_path)

    sampler = Sampler(data_path, sheet_name=sheet_name)
    sampler.write_sampler(model_path, steps=steps)

    make_figs(data_path=data_path, folder=folder_path, model=model_path, total_steps=steps)


if __name__ == '__main__':
    main()
