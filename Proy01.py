""".
Nombre: Proyecto 1

Descripción:

Autores:

    Gonzalez, Pablo 13-10575
    Rodríguez, Mariano 12-10892

Última modificación:    20/05/2016
"""

# -------------------------- LIBRERIAS ---------------------------------------#
import os


# ------------------------- INICIO CLASES ------------------------------------#


class ADNDoble():
    """ Cadena de codones dobles."""
    def __init__(self, cadena):
        """ Define la cadena de la clase y el largo de la misma."""
        self.largo = len(cadena)
        self.cadena = self.complementar(cadena)

    def complementar(self, arreglo):
        """ Complementa una cadena simple, para crear una cadena doble."""
        for i in range(len(arreglo)):
            if arreglo[i] == 'T' or arreglo[i] == 'TA':
                arreglo[i] = 'TA'
            elif arreglo[i] == 'A' or arreglo[i] == 'AT':
                arreglo[i] = 'AT'
            elif arreglo[i] == 'C' or arreglo[i] == 'CG':
                arreglo[i] = 'CG'
            elif arreglo[i] == 'G' or arreglo[i] == 'GC':
                arreglo[i] = 'GC'
            else:
                print("Codón '{}' es inválido, se va a extraer de la cadena."
                      .format(arreglo[i]))
                del arreglo[i]

        return arreglo

    def zip(self, cadenasimple):
        """ Recibe una cadena simple, la complementa y crea la cadena doble."""
        # Copia el contenido de cadenasimple en arreglo
        arreglo = cadenasimple[:]
        # Complementa el arreglo
        self.complementar(arreglo)
        # Guarda el nuevo arreglo en la clase
        self.cadena = arreglo

    def unzip(self, elemento):
        """ Separa una cadena doble en dos cadenas simples."""

        # Inicializa variables
        arreglo1 = ['' for x in range(len(elemento))]
        arreglo2 = ['' for x in range(len(elemento))]

        # Separa la cadena doble en dos simples
        for i in range(len(elemento)):
            arreglo1[i] = elemento[i][0]
            arreglo2[i] = elemento[i][1]

        return arreglo1, arreglo2

    def mitosis(self):
        """ Crea una copia de la doble cadena mediante mitosis."""
        arreglo1, arreglo2 = self.unzip()
        # Complementa las dos cadenas nuevas
        self.mitosis1 = self.complementar(arreglo1)
        self.mitosis2 = self.complementar(arreglo2)

    def buscar(self, cadena, sub_cadena):
        """.
        Recibe una sub_cadena y busca tanto dicha sub_cadena como su
        complemento en ambas cadenas de un ADNDoble."""

        # [subcadena, complemento_subcadena]
        si_cadena = [-1, -1]
        si_complemento = [-1, -1]

        # Complementa la sub_cadena
        sub_cadena = self.complementar(sub_cadena)
        print(sub_cadena)

        # Ciclo de búsqueda en el ADN Doble
        for k in range(2):  # Busca tanto la sub_cadena como su complemento
            for i in range(len(cadena)):
                if cadena[i][0] == sub_cadena[0][k]:
                    for j in range(1, len(sub_cadena)):
                        if cadena[i+j][0] == sub_cadena[j][k] and j == len(sub_cadena) -1:
                            si_cadena[k] = i
                            break
                        elif cadena[i+j][0] == sub_cadena[j][k] and j < len(sub_cadena) -1:
                            continue
                        elif cadena[i+j][0] != sub_cadena[j][k]:
                            si_cadena[k] = -1
                # Si consigue la primera ocurrencia rompe el ciclo
                elif si_cadena[k] != -1:
                    break

        # Ciclo de búsqueda en el complemento del ADN Doble
        for k in range(2):  # Busca tanto la sub_cadena como su complemento
            for i in range(len(cadena)):
                if cadena[i][1] == sub_cadena[0][k]:
                    for j in range(len(sub_cadena)):
                        if cadena[i+j][1] == sub_cadena[j][k] and j == len(sub_cadena) -1:
                            si_complemento[k] = i
                            break
                        elif cadena[i+j][1] == sub_cadena[j][k]:
                            continue
                        elif cadena[i+j][1] != sub_cadena[j][k]:
                            si_complemento[k] = -1
                # Si consigue la primera ocurrencia rompe el ciclo
                elif si_cadena[k] != -1:
                    break

        # Imprime los resultados de la búsqueda
        self.imprimir(cadena)
        print(chr(27) + "[1;37m" + "Se encontraron:" + chr(27) + "[0m")
        print(" La subcadena en la posición: {} de la cadena "
              "y/o {} del complemento"
              .format(si_cadena[0], si_complemento[0]))
        print(" El complemento de la subcadena en la posición: " +
              "{} de la cadena y/o {} del complemento"
              .format(si_cadena[1], si_complemento[1]))

    def imprimir(self, elemento):
        """Escribe en forma de pares."""
        # Imprime en pares en una sola línea.
        for i in elemento:
            print("(" + (chr(27)+"[0;34m" + "{}".format(i[0])) + (chr(27) +
                  "[0;31m" + "{}".format(i[1])) + (chr(27) + "[0m" + ")"),
                  end="")
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
        self.arnt = ARNt(arn)
        return self.arnt


class ARNt():
    """ Cadena de codones triples."""

    proteinas = {  # Diccionario Gxx
                 'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
                 'GCU': 'Val', 'GCC': 'Val', 'GCA': 'Val', 'GCG': 'Val',
                 'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
                 'GGU': 'Gly', 'GGC': 'Gly', 'GGG': 'Gly', 'GGA': 'Gly',
                   # Diccionario Uxx
                 'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
                 'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
                 'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
                 'UGU': 'Cys', 'UGC': 'Cys', 'UGG': 'Stop', 'UGA': 'Trp',
                   # Diccionario Cxx
                 'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
                 'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
                 'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
                 'CGU': 'Arg', 'CGC': 'Arg', 'CGG': 'Arg', 'CGA': 'Arg',
                   # Diccionario Axx
                 'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
                 'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
                 'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
                 'AGU': 'Arg', 'AGC': 'Arg', 'AGG': 'Ser', 'AGA': 'Ser'}

    def __init__(self, cadena):
        """ Define la cadena de la clase y el largo de la misma."""
        self.largo = len(cadena)
        self.cadena = cadena

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

        # Arreglo que contiene la cadena principal en trios
        prot = []
        # Rango de acción del for
        rango = (self.largo//3) * 3
        # Une la cadena de tres en tres
        for i in range(0, rango, 3):
            temp = self.cadena[i] + self.cadena[i+1] + self.cadena[i+2]
            prot.append(temp)

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
                        if j == len(prot):
                            print("basura")
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

        self.proteina = Proteina(traducidos)
        return self.proteina


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
    leidos = []  # Crea el arreglo que contiene las cadenas simples
    with open(archivo, "r") as fd:   # Se abre el archivo para lectura
        for line in fd:  # Se lee linea por linea
            line = line.rstrip()  # Se elimina el salto-de-linea al final de la linea
            line = line.split("\t")
            leidos.append(line)  # llena el arreglo de las cadenas simples

    #  Cada cadena simple de leidos, la complementa y crea la cadena doble
    for i in range(len(leidos)):
        chains = ADNDoble(leidos[i])
        chains.zip(chains.cadena)
        chains.imprimir()

    return chains


def escribir_archivo(A):
    """ Guarda las cosas en un archivo .txt. ROBADO DEL SUDOKU SIN MODIFICAR"""
    i = 0
    p_guardada = 'A.cadena' + str(i) + '.txt'
    existe = os.path.exists('As/' + p_guardada)
    # Verifica si el archivo a guardar ya existe
    while existe is True:
        i += 1
        p_guardada = 'A' + str(i) + '.txt'
        existe = os.path.exists('As/' + p_guardada)
    # Guarda las variables en el archivo (ARREGLAR)
    with open('As/' + p_guardada, 'w') as archivo:
        archivo.write(str(A.modo))
        archivo.write('\n')
        for i in range(A.modo):
            for j in range(A.modo):
                archivo.write(A.tablero_juego[i][j])
            archivo.write('\n')
        archivo.write(A.jugador)
        archivo.write('\n')
        archivo.write(str(A.nivel))
        archivo.write('\n')
        archivo.write(str(A.niveles))
        archivo.write('\n')
        archivo.write(str(A.puntaje))
        archivo.write('\n')
        archivo.write(str(A.tiempo[0]))
        archivo.write('\n')
        archivo.write(str(A.tiempo[1]))
        archivo.write('\n')
        archivo.write(str(A.errores))
        archivo.write('\n')
        archivo.write(str(A.nayuda))
        archivo.write('\n')
        archivo.write(str(A.jugadas))
        archivo.write('\n')
        archivo.write(A.ruta)
        archivo.write('\n')
        archivo.write(A.archivo)
        archivo.write('\n')
        archivo.closed


# ------------------------------ FIN FUNCIONES -------------------------------#

# ------------------------------ PRINCIPAL -----------------------------------#
lista = ['T', 'A', 'C', 'A', 'A', 'A', 'A', 'G', 'A', 'A', 'T', 'A', 'A', 'C',
         'A', 'A', 'T', 'C', 'T', 'A', 'C', 'T', 'A', 'A', 'C', 'T', 'T', 'C',
         'C', 'T', 'T', 'A', 'C', 'A', 'C', 'C']
simple = ADNSimple(lista)
simple.transcribir()
simple.arnt.traducir()

doble = ADNDoble(lista)
doble.buscar(doble.cadena, ['A', 'A', 'A'])


# lectura_cadena('datos2.txt')
