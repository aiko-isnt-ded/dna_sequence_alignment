# Imports
import numpy as np
import time
import matplotlib.pyplot as plt
import seaborn as sns

### NEEDLEMAN-WUNSCH

def needleman_wunsch(seq1, seq2, MATCH_SCORE, MISMATCH_PENALTY, GAP_PENALTY):

    start = time.perf_counter()

    def similarity(a,b, MATCH_SCORE, MISMATCH_PENALTY):
        if a == b:
            return MATCH_SCORE
        else:
            return MISMATCH_PENALTY

    # Creating the empty matrix
    nwm = np.zeros((len(seq1) + 1,len(seq2) + 1))

    # Assigning starter values
    for i in range(nwm.shape[0]):
        for j in range(nwm.shape[1]):
            nwm[i,0] = GAP_PENALTY * i
            nwm[0,j] = GAP_PENALTY * j   

    #Filling out the rest of the matrix
    for i in range(1,nwm.shape[0]):
        for j in range(1, nwm.shape[1]):            
            top_value = nwm[i-1,j] + GAP_PENALTY
            left_value = nwm[i, j-1] + GAP_PENALTY

            if seq1[i-1] == seq2[j-1]:
                diagonal_value = nwm[i-1, j-1] + MATCH_SCORE
            else:
                diagonal_value = nwm[i-1, j-1] + MISMATCH_PENALTY

            nwm[i,j] = max(top_value, left_value, diagonal_value)
    
    #Traceback

    sequence1 = ''
    sequence2 = ''
    path_x = []
    path_y = []

    i, j = len(seq1), len(seq2)

    while (i > 0 or j > 0):
        path_x.append(i)
        path_y.append(j)

        if (i > 0 and j > 0 and nwm[i,j] == nwm[i-1,j-1] + similarity(seq1[i-1], seq2[j-1], MATCH_SCORE, MISMATCH_PENALTY)):
            sequence1 = seq1[i-1] + sequence1
            sequence2 = seq2[j-1] + sequence2
            i -= 1
            j -= 1
        elif (i > 0 and nwm[i,j] == nwm[i-1, j] + GAP_PENALTY):
            sequence1 = seq1[i-1] + sequence1
            sequence2 = '-' + sequence2
            i -= 1
        else:
            sequence1 = '-' + sequence1
            sequence2 = seq2[j-1] + sequence2
            j -= 1
    
    score = nwm[len(seq1), len(seq2)]
    time_taken = time.perf_counter() - start # The time taken to visualize the graph is not considered as time to complete the algorithm, so the timer stops here

    plt.figure(figsize=(8,8))

    sns.heatmap(nwm, xticklabels=[], yticklabels=[])

    plt.plot([x + 0.5 for x in reversed(path_x)],
             [y + 0.5 for y in reversed(path_y)],
             linewidth = 2, color='green')
    
    plt.title("Needleman-Wunsch")
    plt.show()

    return nwm, sequence2, sequence1, score, time_taken

### SMITH-WATERMAN
def smith_waterman(seq1, seq2, match_score, mismatch_score, gap_penalty):
    # Iniciar cronómetro
    Inicio = time.perf_counter() 

    # Longitudes de las secuencias
    n, m = len(seq1), len(seq2)

    # Matriz de puntuaciones (n+1)x(m+1) inicializada en 0
    score_matrix = np.zeros((n+1, m+1), dtype=int)

    # Guardar máxima puntuación y su posición
    max_score = 0
    max_i = max_j = 0

    # Llenar la matriz
    for i in range(1, n+1):
        for j in range(1, m+1):
            # Puntaje match/mismatch
            s = match_score if seq1[i-1] == seq2[j-1] else mismatch_score

            match = score_matrix[i-1, j-1] + s
            delete = score_matrix[i-1, j] + gap_penalty
            insert = score_matrix[i, j-1] + gap_penalty

            # Tomar el máximo entre 0 y los tres movimientos
            score_matrix[i, j] = max(0, match, delete, insert)

            # Guardar la mejor celda
            if score_matrix[i, j] > max_score:
                max_score = score_matrix[i, j]
                max_i, max_j = i, j

    # Traceback
    seq_ord, seq_ord2 = [], []
    path_x, path_y = [], []
    i, j = max_i, max_j

    while score_matrix[i, j] > 0:
        path_x.append(i)
        path_y.append(j)

        s = match_score if seq1[i-1] == seq2[j-1] else mismatch_score

        if score_matrix[i, j] == score_matrix[i-1, j-1] + s:
            seq_ord.append(seq1[i-1])
            seq_ord2.append(seq2[j-1])
            i -= 1; j -= 1
        elif score_matrix[i, j] == score_matrix[i-1, j] + gap_penalty:
            seq_ord.append(seq1[i-1])
            seq_ord2.append('-')
            i -= 1
        else:
            seq_ord.append('-')
            seq_ord2.append(seq2[j-1])
            j -= 1

    seq_ord = ''.join(reversed(seq_ord))
    seq_ord2 = ''.join(reversed(seq_ord2))
    Fin = time.perf_counter()

    # Imprimir los resultados
    print("First Sequence Aligned:", seq_ord)
    print("Second Sequence Aligned:", seq_ord2)
    print("Alignment Score:", max_score)
    print(f"Time Taken: {Fin - Inicio:.4f} seconds")

    # Matriz de confusión con el traceback
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(score_matrix, fmt='d', cmap='YlGnBu', xticklabels=['-']+list(seq2), yticklabels=['-']+list(seq1), ax=ax)
    ax.plot([y+0.5 for y in reversed(path_y)],
            [x+0.5 for x in reversed(path_x)],
            lw=2, color='blue')
    
    ax.scatter(max_j+0.5, max_i+0.5, color='red', s=50, label='Inicio')
    ax.set_title("Matriz Smith-Waterman (alineamiento local)")
    ax.set_ylabel("Secuencia 1")
    ax.set_xlabel("Secuencia 2")
    ax.legend(loc='upper right')
    plt.show()

    return seq_ord, seq_ord2, max_score, score_matrix # Esto es solo si se quiere almacenar el valor en algún lado :)


### Execution

if __name__ == "__main__":
    print('This program allows you to use the Needleman-Wunsch and Smith-Waterman algorithms')
    algo_choice = int(input("Do you want to use Needleman-Wunsch[0] or Smith-Waterman[1]? "))
    first_sequence = str(input("Please copy your first sequence: "))
    second_sequence = str(input("Please copy your second sequence: "))
    user_match = int(input("What is your match score? "))
    user_mismatch = int(input("What is your mismatch penalty? "))
    user_gap = int(input("What is your gap penalty? "))
    print()

    if algo_choice == 0:
        algo_results = needleman_wunsch(first_sequence, second_sequence, user_match, user_mismatch, user_gap)
        print('Here are your results')
        print(f"First Sequence: {algo_results[1]}")
        print(f"Second Sequence: {algo_results[2]}")
        print(f'Alignment Score: {algo_results[3]}')
        print(f"Time Taken: {algo_results[4]:.10f}")
    
    elif algo_choice == 1:
        smith_waterman(first_sequence, second_sequence, user_match, user_mismatch, user_gap)

    else:
        print('Unable to process your entry. Try again.')