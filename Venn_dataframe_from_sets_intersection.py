class VennIntersectionsEnhanced:
    def __init__(self, **kwargs):
        """
        Initialize the class with named sets.
        Example usage:
            VennIntersectionsEnhanced(set1=setA, set2=setB, set3=setC)
        """
        self.sets = kwargs
        self.intersections = {}
        self.calculate_intersections()

    def calculate_intersections(self):
        # Reset intersections
        self.intersections = {}
        keys = list(self.sets.keys())
        n = len(keys)

        # Iterate through all possible combinations of the sets
        for i in range(1, n + 1):
            for combo in combinations(keys, i):
                intersected_set = set.intersection(*[self.sets[key] for key in combo])
                if intersected_set:
                    # Create a key representing the combination
                    key = '_'.join(combo)
                    self.intersections[key] = intersected_set

    def get_intersections(self):
        return self.intersections

    def print_intersections(self):
        for key, value in self.intersections.items():
            print(f"{key}: {value}")

    def intersections_to_dataframe(self):
        # Convert intersections to DataFrame
        # Find the longest list
        max_length = max(len(v) for v in self.intersections.values())

        # Pad shorter lists with NaN
        padded_intersections = {k: list(v) + [None]*(max_length - len(v)) for k, v in self.intersections.items()}

        # Create DataFrame
        return pd.DataFrame.from_dict(padded_intersections, orient='index').transpose()

    
    
# Example usage with the enhanced class
# Example usage
setA = {'gene1', 'gene2', 'gene3'}
setB = {'gene2', 'gene3', 'gene4'}
setC = {'gene3', 'gene4', 'gene5'}
setD = {'gene6', 'gene7'}

venn_enhanced = VennIntersectionsEnhanced(group1=setA, group2=setB, group3=setC, group4=setD)

# Print the intersections
venn_enhanced.print_intersections()

# Convert intersections to DataFrame and display it
df_enhanced = venn_enhanced.intersections_to_dataframe()
df_enhanced

