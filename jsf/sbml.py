import libsbml
import math


def read_sbml(sbml_xml):
    reader = libsbml.SBMLReader()
    doc = reader.readSBML(sbml_xml)
    model = doc.getModel()
    num_species = model.getNumSpecies()
    num_reactions = model.getNumReactions()

    assert model.getCompartment(0).getId() == "main", "Main compartment not found"

    if num_species == 1:
        x0 = model.getSpecies(0).getInitialAmount()
        stoichiometry = []
        rate_param = []
        reaction_details = []
        for r_ix in range(num_reactions):
            reaction = model.getReaction(r_ix)
            stoichiometry.append(reaction.getReactant(0).getStoichiometry())
            rate_param.append(reaction.getKineticLaw().getParameter(0).getValue())
            reactants = [
                (reactant.getSpecies(), reactant.getStoichiometry())
                for reactant in reaction.getListOfReactants()
            ]
            products = [
                (product.getSpecies(), product.getStoichiometry())
                for product in reaction.getListOfProducts()
            ]
            reaction_details.append(
                {
                    "reaction_id": reaction.getId(),
                    "reactants": reactants,
                    "products": products,
                    "rate_parameter": reaction.getKineticLaw()
                    .getParameter(0)
                    .getValue(),
                }
            )

        def rates(x, t):
            return [
                rate_param[i] * x[0] ** stoichiometry[i] for i in range(num_reactions)
            ]

        def _zero_if_missing(x, y):
            tmp = x[y]
            return 0 if tmp == [] else tmp[0][-1]

        nu_reactants = [[_zero_if_missing(r, "reactants")] for r in reaction_details]
        nu_products = [[_zero_if_missing(r, "products")] for r in reaction_details]

        # # In the nu-matrix each row is a reaction and each column
        # # describes the number items of that species used in the
        # # reaction.
        stoich = {
            "nu": [
                [a - b for a, b in zip(r1, r2)]
                for r1, r2 in zip(nu_products, nu_reactants)
            ],
            "DoDisc": [1],
            "nuReactant": nu_reactants,
            "nuProduct": nu_products,
        }
    elif num_species > 1:
        x0 = [model.getSpecies(i).getInitialAmount() for i in range(num_species)]
        r_ix_map = {model.getSpecies(i).getId(): i for i in range(num_species)}

        nu_reactants = [[0] * num_species for r in range(num_reactions)]
        nu_products = [[0] * num_species for r in range(num_reactions)]

        stoichiometry = []
        rate_param = []
        reaction_details = []
        for r_ix in range(num_reactions):
            reaction = model.getReaction(r_ix)
            rate_param.append(reaction.getKineticLaw().getParameter(0).getValue())
            reactants = [
                (reactant.getSpecies(), reactant.getStoichiometry())
                for reactant in reaction.getListOfReactants()
            ]
            products = [
                (product.getSpecies(), product.getStoichiometry())
                for product in reaction.getListOfProducts()
            ]
            reaction_details.append(
                {
                    "reaction_id": reaction.getId(),
                    "reactants": reactants,
                    "products": products,
                    "rate_parameter": reaction.getKineticLaw()
                    .getParameter(0)
                    .getValue(),
                }
            )

            for r in reaction.getListOfReactants():
                species_id = r.getSpecies()
                species_ix = r_ix_map[species_id]
                nu_reactants[r_ix][species_ix] = r.getStoichiometry()

            for p in reaction.getListOfProducts():
                species_id = p.getSpecies()
                species_ix = r_ix_map[species_id]
                nu_products[r_ix][species_ix] = p.getStoichiometry()

        def rates(x, t):
            return [
                rate_param[i]
                * math.prod(x[j] ** nu_reactants[i][j] for j in range(num_species))
                for i in range(num_reactions)
            ]

        # # In the nu-matrix each row is a reaction and each column
        # # describes the number items of that species used in the
        # # reaction.
        stoich = {
            "nu": [
                [a - b for a, b in zip(r1, r2)]
                for r1, r2 in zip(nu_products, nu_reactants)
            ],
            "DoDisc": [1 for _ in range(num_species)],
            "nuReactant": nu_reactants,
            "nuProduct": nu_products,
        }
    else:
        raise ValueError("No species found in the model")

    return x0, rates, stoich
