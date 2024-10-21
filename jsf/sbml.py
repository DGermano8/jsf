import libsbml

def read_sbml(sbml_xml):
    reader = libsbml.SBMLReader()
    doc = reader.readSBML(sbml_xml)
    model = doc.getModel()

    assert model.getCompartment(0).getId() == 'main', 'Main compartment not found'

    if model.getNumSpecies() == 1:
        x0 = model.getSpecies(0).getInitialAmount()
        num_reactions = model.getNumReactions()
        stoichiometry = []
        rate_param = []
        reaction_details = []
        for r_ix in range(num_reactions):
            reaction = model.getReaction(r_ix)
            stoichiometry.append(reaction.getReactant(0).getStoichiometry())
            rate_param.append(reaction.getKineticLaw().getParameter(0).getValue())
            reactants = [(reactant.getSpecies(), reactant.getStoichiometry()) for reactant in reaction.getListOfReactants()]
            products = [(product.getSpecies(), product.getStoichiometry()) for product in reaction.getListOfProducts()]
            reaction_details.append({
                "reaction_id": reaction.getId(),
                "reactants": reactants,
                "products": products,
                "rate_parameter": reaction.getKineticLaw().getParameter(0).getValue(),
            })

        def rates(x, t):
                return [rate_param[i] * x[0]**stoichiometry[i] for i in range(num_reactions)]

        def _zero_if_missing(x, y):
            tmp = x[y]
            return 0 if tmp == [] else tmp[0][-1]

        nu_reactants = [[_zero_if_missing(r,'reactants')] for r in reaction_details]
        nu_products = [[_zero_if_missing(r,'products')] for r in reaction_details]

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
    else:
        raise ValueError('Multiple species not implemented yet')


    return x0, rates, stoich
