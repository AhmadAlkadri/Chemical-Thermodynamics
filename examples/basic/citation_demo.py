"""Demonstration of citation retrieval for thermodynamic data."""

from chemthermo import Component, cite

def main():
    print("=== ChemThermo Citation Demo ===\n")

    # 1. Direct citation lookup
    print("--- Direct Lookup ---")
    print(f"Citation for Methane Tc: {cite('Methane', 'Tc')}")
    print(f"Citation for Water MW:   {cite('Water', 'MW')}")
    print()

    # 2. Component object access
    print("--- Component Object Access ---")
    eth = Component.from_database("Ethanol")
    print(f"Component: {eth.name}")
    print(f"Formula:   {eth.formula}")
    print(f"Tc Source: {eth.get_citation('Tc')}")
    
    if eth.antoine:
        print(f"Antoine:   {eth.get_citation('antoine')}")
    else:
        print("Antoine:   (Not available)")
        
    print()

    # 3. Handling missing citations
    print("--- Missing/Unavailable ---")
    # This might show "Legacy ..." if everything is legacy, or "No citation available"
    print(f"Citation for non-existent property 'foobar': {eth.get_citation('foobar')}")

if __name__ == "__main__":
    main()
