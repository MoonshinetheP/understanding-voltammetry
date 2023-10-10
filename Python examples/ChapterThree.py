import numpy as np

class ChThree:

    def __init__(self, oEini = 20, oEfin = -20, osr = 100, odX = 4E-4, odE = 0.02):
        
        self.oEini = oEini
        self.oEfin = oEfin
        self.osr = osr
        self.odX = odX
        self.odE = 0.02

        self.odT = self.odE / self.osr
        self.omaxT = 2 * np.abs(self.oEfin - self.oEini) / self.osr
        self.omaxX = 6 * np.sqrt(self.omaxT)
        self.n = int(1 + self.omaxX / self.odX)
        self.m = int(self.omaxT / self.odT)

        self.l = self.odT / (self.odX * self.odX)
        self.a = -self.l
        self.b = 2 * self.l + 1
        self.g = -self.l

        self.gmod = np.zeros((self.n - 1,))
        self.delta = np.ones((self.n - 1,))
        self.dmod = np.zeros((self.n - 1,))
        self.C = np.ones((self.n))

        self.gmod[0] = 0
        for i in range(1, self.n -1, 1): # same
            self.gmod[i] = self.g / (self.b - self.gmod[i-1] * self.a) 

        self.oE = self.oEini

        self.E = np.array([])
        self.flux = np.array([])

        for k in range(0, self.m, 1):
            if k < self.m/2:
                self.oE -= self.odE
            else:
                self.oE += self.odE
        
            self.dmod[0] = 1 / (1 + np.exp(-self.oE))
            for x in range(1, self.n - 1):
                self.dmod[x] = (self.delta[x] - self.dmod[x - 1] * self.a) / (self.b - self.gmod[x - 1] * self.a)

            self.C[self.n - 1] = 1
            for y in range(self.n - 2, -1, -1):
                self.C[y] = self.dmod[y] - self.gmod[y] * self.C[y + 1]

                self.delta[y] = self.C[y]

            
            self.E = np.append(self.E, self.oE)
            self.flux = np.append(self.flux, (-(-self.C[2] + 4 * self.C[1] - 3 * self.C[0])) / (2 * self.odX))
        #self.output = zip(self.E, self.flux)

    def results(self):
        return self.flux


if __name__ == '__main__':
    instance = ChThree()
    data = 'C:/Users/SLinf/Documents/data.txt'
    with open(data, 'w') as file:
        for ix in instance.flux:
            file.write(str(ix) + '\n')