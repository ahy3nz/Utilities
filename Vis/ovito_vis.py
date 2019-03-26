import ovito


def render(painter, **args):
    node = ovito.dataset.selected_node
    node.compute()
    atypes = node.source.particle_properties.particle_type.type_list
    node.output.bonds.display.width = 0.3
    for atype in atypes:
        atype.radius=0.1
    print("DIDSTUFF")
