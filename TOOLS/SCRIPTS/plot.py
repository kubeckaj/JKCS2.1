import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Plot data from coloumn divided text file.')
    parser.add_argument('filename', type=str, help='Path to the text file containing the data.')
    parser.add_argument('-gui', action='store_true', help='Show matplotlib GUI')
    parser.add_argument('-xscale', type=str, default='linear', choices=['linear', 'log', 'symlog', 'logit'], help='Scale type for the X axis.')
    parser.add_argument('-yscale', type=str, default='linear', choices=['linear', 'log', 'symlog', 'logit'], help='Scale type for the Y axis.')
    parser.add_argument('-xname', type=str, default='X', help='Set x-axis name')
    parser.add_argument('-yname', type=str, default='Y', help='Set y-axis name')
    parser.add_argument('-title', type=str, default='Data Plot', help='Title of the plot.')
    parser.add_argument('-plot_type', type=str, default='line', choices=['line', 'scatter', 'bar', 'hist'], help='Type of plot.')
    parser.add_argument('-linestyle', type=str, default='-', choices=[':', '--', '-.'], help='Linestyle')
    parser.add_argument('-regression', action='store_true', default=False, help='Plot linear regression line.')
    parser.add_argument('-grid', action='store_true', help='Show grid in plot')
    parser.add_argument('-x', type=int, default=1, help='Column number for X axis data.')
    parser.add_argument('-y', type=int, default=2, help='Column number for Y axis data.')
    return parser.parse_args()

def load_data(filename, x_col=None, y_col=None):
    usecols = (x_col - 1, y_col - 1) if x_col and y_col else None
    data = np.loadtxt(filename, usecols=usecols)
    return data

def plot_data(data, args):
    plt.figure()
    if args.plot_type == 'line':
        plt.plot(data[:, 0], data[:, 1], linestyle=args.linestyle)
    elif args.plot_type == 'scatter':
        plt.scatter(data[:, 0], data[:, 1], label='Data')
        if args.regression:
            fit_coefficients = np.polyfit(data[:, 0], data[:, 1], 1)
            fit_function = np.poly1d(fit_coefficients)
            x_values_for_fit = np.linspace(min(data[:, 0]), max(data[:, 0]), num=100)  # Change 0 to 1 if you want to spah over y range
            plt.plot(x_values_for_fit, fit_function(x_values_for_fit), color='red', linestyle=args.linestyle, label=f'Linear fit: y={fit_coefficients[0]:.2f}x+{fit_coefficients[1]:.2f}')
            plt.legend()
    elif args.plot_type == 'bar':
        plt.bar(data[:, 0], data[:, 1])
    elif args.plot_type == 'hist':
        plt.hist(data)

    plt.xscale(args.xscale)
    plt.yscale(args.yscale)
    plt.title(args.title)
    plt.xlabel(f'{args.xname}')
    plt.ylabel(f'{args.yname}')
    plt.grid(args.grid)
    if args.gui:
        plt.show()
    else:
        plt.savefig(f"plot_{args.xname}_{args.yname}.svg", dpi=600)

def main():
    args = parse_args()

    try:
        # Depending on the plot type, different data loading might be needed
        if args.plot_type in ['line', 'scatter', 'bar']:
            data = load_data(args.filename, args.x, args.y)
        elif args.plot_type == 'hist':
            data = load_data(args.filename)
        plot_data(data, args)
    except Exception as e:
        print(f"Error: {str(e)}")
        exit()

if __name__ == '__main__':
    main()

