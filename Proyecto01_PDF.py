""".
Nombre: Proyecto 1

Descripción:

Autores:

    Gonzalez, Pablo 13-10575
    Rodríguez, Mariano 12-10892

Última modificación:    20/05/2016
"""

# Librerias
import os


# ------------------------- INICIO CLASES ------------------------------------#

class Contenedora():
    """.
    Objeto que cotiene todos los datos que se generan o usan en el programa."""
    ADNBasura = []
    ARNBasura = []
    proteina = []

    def __init__(self, dato_simple):
        # Crea el ADNSimple
        self.simple = ADNSimple(dato_simple)
        print('Cadena simple: ', self.simple.cadena)
        # Crea el ADNDoble
        self.doble = ADNDoble(self.simple.cadena)
        self.doble.cadena = self.doble.zip()
        # Crea el ARNt
        self.arnt = ARNt(self.simple.transcribir())

    def escribir(self, archivo=None):
        """ Crea los archivos de texto con los resultados de las corridas."""
        global ADNBasura, ARNBasura
        if archivo is None:
            i = 0
            s_guardada = 'resultados/organismo' + str(i) + '.txt'
            existe = os.path.exists(s_guardada)
            print(existe)
            while existe is True:
                i += 1
                s_guardada = 'resultados/organismo' + str(i) + '.txt'
                existe = os.path.exists(s_guardada)
            print(s_guardada)

            with open(s_guardada, 'w') as f:
                """for i in range(len(self.cadena)):
                    f.write(self.cadena[i])
                    f.write('\t')
                f.write('\n')"""
                f.write("ADNSimple:  " + str(self.simple.cadena))
                f.write('\n')
                f.write("ADNDoble:   " + str(self.doble.cadena))
                f.write('\n')
                f.write("ARNt:   " + str(self.arnt.cadena))
                f.write('\n')
                f.write("Proteinas:   " + str(self.proteina.proteinas))
                f.write('\n')
                f.write("ADNBasura:   " + str(ADNBasura))
                f.write('\n')
                f.write("ARNBasura:   " + str(ARNBasura))
                f.write('\n')
        else:
            with open(archivo + '.txt', 'a') as f:
                """for i in range(len(self.cadena)):
                    f.write(self.cadena[i])
                    f.write('\t')
                f.write('\n')"""
                f.write(str(self.simple.cadena))
                f.write('\n')
                f.write(str(self.doble.cadena))
                f.write('\n')
                f.write(str(self.arnt.cadena))
                f.write('\n')
                f.write(str(self.proteina.proteinas))
                f.write('\n')
                f.write(str(ADNBasura))
                f.write('\n')
                f.write(str(ARNBasura))
                f.write('\n')


class ADNDoble():
    """ Cadena de codones dobles."""
    def __init__(self, cadena):
        """ Define la cadena de la clase y el largo de la misma."""
        self.largo = len(cadena)
        self.cadena = cadena

    def complementar(self, arreglo):
        """ Complementa una cadena simple, para crear una cadena doble."""
        for i in range(len(arreglo)):
            if arreglo[i] == 'T':
                arreglo[i] = 'TA'
            elif arreglo[i] == 'A':
                arreglo[i] = 'AT'
            elif arreglo[i] == 'C':
                arreglo[i] = 'CG'
            elif arreglo[i] == 'G':
                arreglo[i] = 'GC'
            else:
                print("Codón '{}' es inválido, se va a extraer de la cadena."
                      .format(arreglo[i]))
                del arreglo[i]

        return arreglo

    def zip(self):
        """ Recibe una cadena simple, la complementa y crea la cadena doble."""
        # Copia el contenido de cadenasimple en arreglo
        arreglo = self.cadena[:]
        # Complementa el arreglo
        self.complementar(arreglo)
        # Guarda el nuevo arreglo en la clase
        self.cadena = arreglo

        return arreglo

    def unzip(self):
        """ Separa una cadena doble en dos cadenas simples."""

        # Inicializa variables
        arreglo1 = ['' for x in range(self.largo)]
        arreglo2 = ['' for x in range(self.largo)]

        # Separa la cadena doble en dos simples
        for i in range(self.largo):
            arreglo1[i] = self.cadena[i][0]
            arreglo2[i] = self.cadena[i][1]

        return arreglo1, arreglo2

    def mitosis(self):
        """NP Crea una copia de la doble cadena mediante mitosis."""
        arreglo1, arreglo2 = self.unzip()
        # Complementa las dos cadenas nuevas
        self.mitosis1 = self.complementar(arreglo1)
        self.mitosis2 = self.complementar(arreglo2)

    def buscar(self, cadena, sub_cadena):
        """Recibe una sub_cadena y busca tanto dicha sub_cadena como su
        complemento en ambas cadenas de un ADNDoble."""

        # [subcadena, complemento_subcadena]
        si_cadena = [-1, -1]
        si_complemento = [-1, -1]

        # Complementa la sub_cadena
        sub_cadena = self.complementar(sub_cadena)
        print(sub_cadena)
        # print(cadena)

        # Ciclo de búsqueda en el ADN Doble
        for k in range(2):  # Busca tanto la sub_cadena como su complemento
            for i in range(len(cadena)):
                if cadena[i][0] == sub_cadena[0][k]:
                    for j in range(1, len(sub_cadena)):
                        if cadena[i+j][0] == sub_cadena[j][k] and (j == len(sub_cadena) -1):
                            si_cadena[k] = i
                            break
                        elif cadena[i+j][0] == sub_cadena[j][k] and j < len(sub_cadena) -1:
                            pass
                        elif cadena[i+j][0] != sub_cadena[j][k]:
                            si_cadena[k] = -1
                            break
                # Si consigue la primera ocurrencia rompe el ciclo
                if si_cadena[k] != -1:
                    break

        # Ciclo de búsqueda en el complemento del ADN Doble
        for k in range(2):  # Busca tanto la sub_cadena como su complemento
            for i in range(len(cadena)):
                if cadena[i][1] == sub_cadena[0][k]:
                    for j in range(1, len(sub_cadena)):
                        if cadena[i+j][1] == sub_cadena[j][k] and (j == len(sub_cadena) -1):
                            si_complemento[k] = i
                            break
                        elif cadena[i+j][1] == sub_cadena[j][k] and j < len(sub_cadena) -1:
                            pass
                        elif cadena[i+j][1] != sub_cadena[j][k]:
                            si_complemento[k] = -1
                            break
                # Si consigue la primera ocurrencia rompe el ciclo
                if si_complemento[k] != -1:
                    break

        # Imprime los resultados de la búsqueda
        self.imprimir()
        print(chr(27) + "[1;37m" + "Se encontraron:" + chr(27) + "[0m")
        print(" La subcadena en la posición: {} de la cadena "
              "y/o {} del complemento"
              .format(si_cadena[0], si_complemento[0]))
        print(" El complemento de la subcadena en la posición: " +
              "{} de la cadena y/o {} del complemento"
              .format(si_cadena[1], si_complemento[1]))

    def imprimir(self):
        """Escribe en forma de pares."""
        # Imprime en pares en una sola línea.
        for i in self.cadena:
            print("(" + (chr(27)+"[0;34m" + "{}".format(i[0])) + (chr(27) +
                  "[0;31m" + "{}".format(i[1])) + (chr(27) + "[0m" + ")"), end="")
        # Crea la linea al final como separación.
        print('')


# -------------- NP todo hacia abajo


class ADNSimple():
    """ Cadena de codones simples."""

    def __init__(self, cadena):
        """ Define la cadena de la clase y el largo de la misma."""
        self.largo = len(cadena)
        self.cadena = cadena

    def complementar(self, arreglo):
        """ Complementa una cadena simple, para crear una cadena doble."""
        arreglo_comp = ['' for x in range(len(arreglo))]
        for i in range(len(arreglo)):
            if arreglo[i] == 'T':
                arreglo_comp[i] = 'A'
            elif arreglo[i] == 'A':
                arreglo_comp[i] = 'T'
            elif arreglo[i] == 'C':
                arreglo_comp[i] = 'G'
            elif arreglo[i] == 'G':
                arreglo_comp[i] = 'C'
            else:
                print("Codón '{}' es inválido, se va a extraer de la cadena."
                      .format(arreglo[i]))
                del arreglo[i]

        self.complemento = arreglo_comp
        return self.complemento

    def transcribir(self):
        """.
        Devuelve como un ARNt el complemento de la cadena simple,
        pero con la Timina sustituida por Uracilo."""

        # Crea complemento de la cadena original (Futuro ARNt)
        arn = self.complementar(self.cadena)

        # Cambia la Timina por Uracilo
        for i in range(self.largo):
            if arn[i] == 'T':
                arn[i] = 'U'

        # Devuelve el arn
        return arn


class ARNt():
    """ Cadena de codones triples."""

    proteinas = {  # Diccionario Gxx
                 'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
                 'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
                 'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
                 'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
                   # Diccionario Uxx
                 'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
                 'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
                 'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
                 'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp',
                   # Diccionario Cxx
                 'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
                 'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
                 'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
                 'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
                   # Diccionario Axx
                 'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
                 'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
                 'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
                 'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg'}

    def __init__(self, cadena):
        """ Define la cadena de la clase y el largo de la misma."""
        self.largo = len(cadena)
        self.cadena = cadena

    def heap(self):
        """ Crea el arreglo prot a ser ordenado con heap. """

        global ADNBasura, ARNBasura  # Variable global para almacenamiento de trozos basura (en caso de encontrar)

        # Arreglo que contiene la cadena principal en trios
        prot = []
        # Rango de acción del for
        rango = (self.largo//3) * 3
        # Une la cadena de tres en tres
        for i in range(0, rango, 3):
            temp = self.cadena[i] + self.cadena[i+1] + self.cadena[i+2]
            prot.append((temp,0))  # Le agregamos la segunda posición que representa la frecuencia

        prot_heap = prot[:]  # Inicializamos el arreglo a ordenar en heap

        for i in range(len(prot_heap)):  # Contamos las veces que se encuentra el codón en el arreglo
            cont = 0  # Inicio el contador de frecuencia
            for j in range(len(prot_heap)):
                if prot_heap[i][0] == prot_heap[j][0]:
                    cont += 1
                    prot_heap[i] = (prot_heap[i][0],cont)

        ordena_puntero(prot_heap,len(prot_heap))  # Ordena el arreglo con el heapsort (ordena_puntero)

        #  Impresion de los codones y su frecuencia en el arreglo
        pri = ''
        for i in range(len(prot_heap)):
            if pri != prot_heap[i]:
                pri = prot_heap[i]
                print('El ',prot_heap[i][0],' aparece: ',prot_heap[i][1],' veces')

        # Inicializando variables
        x = 0
        traducidos = []
        final = False
        i = 0
        while final is False:
            print('ciclo heap')
            inicio = None
            fin = None
            # Halla si hay una secuencia de inicio en una cadena.
            while i < len(prot_heap):
                # Termina el ciclo si se llega al final del arreglo
                if i == len(prot_heap):
                    final = True
                    break

                if prot_heap[i][0] == 'AUG':
                    inicio = i
                    for j in range(i+1, len(prot_heap) + 1):
                        if j == len(prot_heap):
                            print("BASURA HEAP")
                            break
                        elif self.proteinas[prot_heap[j][0]] == 'Stop':
                            x = j
                            # Traduce en el rango entre inicio y stop
                            for k in range(i+1, j):
                                traducidos.append(self.proteinas[prot_heap[k][0]])
                                # Agrega un '/' al final de la proteina
                                if k == j-1:
                                    traducidos.append('/')
                            i = x
                            break
                    if i == len(prot_heap) - 1:
                        break

                i = i + 1

            final = True

        print('Proteina del heap: ',traducidos)

    def traducir(self):
        """.
        Complementa el ARNt, busca el codón de inicio y devolverá la proteína
        que resulta de encadenar los aminoácidos producto del recorrido de la
        secuencia según la Tabla de aminoácidos, hasta que llega a un stop.
        La traducción comienza con el codón "AUG" que es además de señal de
        inicio (el aminoácido met, metionina). La secuencia de codones
        determina la secuencia de aminoácidos en una proteína en concreto, con
        una estructura y una función específica. Ni el codón de inicio ni el
        codón de stop formarán parte de la secuencia proteica."""

        global ADNBasura, ARNBasura # Variable global para almacenamiento de trozos basura (en caso de encontrar)

        # Todo el proceso del heap
        self.heap()

        # Arreglo que contiene la cadena principal en trios
        prot = []
        # Rango de acción del for
        rango = (self.largo//3) * 3
        # Une la cadena de tres en tres
        for i in range(0, rango, 3):
            temp = self.cadena[i] + self.cadena[i+1] + self.cadena[i+2]
            prot.append(temp)  # Le agregamos la segunda posición que representa la frecuencia

        # Inicializando variables
        x = 0
        traducidos = []
        final = False
        while final is False:
            inicio = None
            fin = None
            # Halla si hay una secuencia de inicio en una cadena.
            for i in range(x, len(prot)+1):
                # Termina el ciclo si se llega al final del arreglo
                if i == len(prot):
                    final = True
                    break

                if prot[i] == 'AUG':
                    inicio = i
                    for j in range(i+1, len(prot)+1):
                        if j == len(prot):  # Encuentra parte basura
                            print("BASURA")
                            ADNBasura.append(A.simple.cadena[3*i:])
                            ARNBasura.append(self.cadena[3*i:])
                        elif self.proteinas[prot[j]] == 'Stop':
                            x = j
                            # Traduce en el rango entre inicio y stop
                            for k in range(i+1, j):
                                traducidos.append(self.proteinas[prot[k]])
                                # Agrega un '/' al final de la proteina
                                if k == j-1:
                                    traducidos.append('/')
                            break

                else:
                    pass

        return traducidos  # ADNBasura, ARNBasura


class Proteina():
    """ Cadena de proteinas, obtenidas de una cadena inicialmente de ARNt."""

    def __init__(self, cadena):
        """ Define la cadena de la clase y el largo de la misma."""
        proteinas = []
        temporal = []
        for i in range(len(cadena)):
            if cadena[i] != '/':
                temporal.append(cadena[i])
            elif cadena[i] == '/':
                # Se agrega la proteina
                proteinas.append(temporal)
                temporal = []

        # Se crean los elementos de la clase
        self.proteinas = proteinas
        self.cantidad_proteinas = len(proteinas)
        self.cadena = cadena


# --------------------------- FIN CLASES -------------------------------------#

# --------------------------- INICIO FUNCIONES -------------------------------#


def lectura_cadena(archivo):
    """ Lee un archivo .txt """
    leidos = []  # Crea el arreglo que contiene las cadenas simples
    with open(archivo, "r") as fd:   # Se abre el archivo para lectura
        for line in fd:                    # Se lee linea por linea
            temp = []
            line = line.rstrip()    # Se elimina el salto-de-linea al final de la linea
            # line = line.split("\t")
            for i in range(len(line)):
                temp.append(line[i])
            leidos.append(temp)  # llena el arreglo de las cadenas simples

    return leidos

    """#  Cada cadena simple  de le complementa y e crea la cadena doble
    for i in range(len(leidos)):
        chains = ADNDoble(leidos[i])
        chains.zip(chains.cadena)
        chains.imprimir()"""


def heapify(A, idx, maxi):

    izquierda = 2 * idx + 1
    derecha = 2 * idx + 2

    # Busca el elemnto mas grande
    if izquierda < maxi and A[izquierda][1] < A[idx][1]:
        pequeño = izquierda
    else:
        pequeño = idx

    if derecha < maxi and A[derecha][1] < A[pequeño][1]:
        pequeño = derecha

    # Algo sobre 'propagate'
    if pequeño != idx:
        A[idx], A[pequeño] = A[pequeño], A[idx]
        heapify(A, pequeño, maxi)


def construye_heap(A, n):
    for i in range(n//2, -1, -1):
        heapify(A, i, n)


def ordena_puntero(A, n):
    construye_heap(A, n)
    for i in range(n-1, 0, -1):
        A[0], A[i] = A[i], A[0]
        heapify(A, 0, i)
    return A


def partition(A, p, r):  # It Worked (con quicksort)
    '''Utiliza un elemento como el pivote.'''
    x = A[r]
    i = p - 1
    for j in range(p, r):
        if len(A[j]) > len(x):
            i += 1
            A[i], A[j] = A[j], A[i]
        elif len(A[j]) == len(x):
            if A[j] <= x:
                i += 1
                A[i], A[j] = A[j], A[i]
    A[i + 1], A[r] = A[r], A[i + 1]
    return i + 1


def quicksort(A, p, r):  # It Worked!
    if p < r:
        q = partition(A, p, r)
        quicksort(A, p, q-1)
        quicksort(A, q+1, r)

    return(A)

# ------------------------------ FIN FUNCIONES -------------------------------#

# ------------------------------ PRINCIPAL -----------------------------------#

# Corrida de simple.txt
datos = lectura_cadena('simple.txt')
for i in range(len(datos)):
    ADNBasura = []
    ARNBasura = []
    A = Contenedora(datos[i])
    A.proteina = A.arnt.traducir()
    A.proteina = Proteina(A.proteina)
    A.escribir()

# Corrida de complejo.txt
datos = lectura_cadena('complejo.txt')
for i in range(len(datos)):
    ADNBasura = []
    ARNBasura = []
    A = Contenedora(datos[i])
    A.proteina = A.arnt.traducir()
    A.proteina = Proteina(A.proteina)
    A.escribir()

# Corrida de aminoacidos.txt
datos = lectura_cadena('aminoacidos.txt')
for i in range(len(datos)):
    ADNBasura = []
    ARNBasura = []
    A = Contenedora(datos[i])
    A.proteina = A.arnt.traducir()
    A.proteina = Proteina(A.proteina)
    A.escribir()
