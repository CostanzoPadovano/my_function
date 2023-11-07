def create_evenn_dataframe(*lists, list_names=None):
    """
    Crea un dataframe per la visualizzazione di un diagramma di Venn.
    
    Parametri:
    *lists: un numero variabile di liste contenenti elementi.
    list_names: una lista di stringhe che rappresentano i nomi delle liste.
    
    Ritorna:
    Un DataFrame con due colonne: 'Element' e 'Group', pronto per un diagramma di Venn.
    """
    
    if list_names is None:
        # Assegna nomi di lista predefiniti se non vengono forniti
        list_names = [f'Lista{i+1}' for i in range(len(lists))]
    elif len(list_names) != len(lists):
        raise ValueError("Il numero dei nomi delle liste deve corrispondere al numero delle liste fornite.")
    
    # Creiamo il dataframe unendo tutte le liste con i rispettivi nomi delle liste
    venn_data = []
    for list_idx, current_list in enumerate(lists):
        # Per ogni elemento nella lista, aggiungi una riga con il nome della lista corrispondente
        venn_data.extend((element, list_names[list_idx]) for element in current_list)
    
    # Creiamo il dataframe
    df_venn = pd.DataFrame(venn_data, columns=['Element', 'Group'])
    
    return df_venn

# Testiamo la funzione con tre liste di esempio
test_list1 = ["A", "B", "C", "D"]
test_list2 = ["B", "C", "E"]
test_list3 = ["C", "D", "F"]

# Creiamo il dataframe
df_test_venn = create_evenn_dataframe(test_list1, test_list2, test_list3, list_names=['Group1', 'Group2', 'Group3'])
df_test_venn.head()  # Mostra le prime righe del dataframe per conferma
