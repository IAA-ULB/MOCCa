import pandas
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    for filename in ('timings','timings2'):
        # clean up the output file:
        lines = []
        with open(filename+'.txt',mode='r') as f:
            for line in f:
                try:
                    i = int(line[0])
                    lines.append(line)
                except ValueError:
                    continue

        with open(filename+'.csv',mode='w') as f:
            for line in lines:
                f.write(line)

        timings = pandas.read_csv('timings.csv',sep=','
                                , names=['rank','np','n','bf','dt'])
        timings = timings.drop('rank', axis=1)
        # print(timings)

        r0 = None
        ir = 0
        sum = pandas.DataFrame()
        for index, row in timings.iterrows():
            if r0 is None:
                r0 = row
            else:
                if  row['np'] == r0['np'] and row['n'] == r0['n'] and row['bf'] == r0['bf']:
                    r0['dt'] += row['dt']
                else:
                    sum[ir] = r0
                    ir += 1
                    r0 = row

        subsets = set()

        sum = sum.transpose()
        for index, row in sum.iterrows():
            # print(row)
            subset = (row['n'],row['bf'])
            subsets.add(subset)

        dataframes = []
        for subset in subsets:
            dataframes.append(sum.loc[(sum['n'] == subset[0]) & (sum['bf'] == subset[1])])
            # print(dataframes[-1])

        for dataframe in dataframes:
            plt.figure()
            x = dataframe['np'].to_numpy()
            idx = np.argsort(x)
            x = x[idx]
            y = dataframe['dt'].to_numpy()[idx]
            plt.plot(x, y, 'o-')
            plt.xlabel('np')
            plt.ylabel('dt')
            n = int(dataframe['n'].iloc[0])
            bf = int(dataframe['bf'].iloc[0])
            title = f"{n}x{n} {bf=}"
            plt.title(title)
            plt.savefig(f"{n}x{n}_{bf}.png")
            dataframe.to_csv(f'{title}.csv', index=False)
            plt.close()