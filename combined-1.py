import bpy
import os
import numpy as np
import sys
import argparse
from astropy import units as u
from astropy import constants as const 

#TO DO: 
#9. subdivide and shade smooth limb sphere, also make it translucent
#10. check if inputfromuser is correct or not
#11. some scripts use eevee, some cycles
#emission in white dwarf
SCALE_FAC_MASS = 5
SCALE_FAC_TEMP = 55
SCALE_FAC_RADIUS = 2

#The following constants have SI units. The output from equations is scaled from SI to units wrt Sun.
M0 = const.M_sun
L0 = const.L_sun
T0 = 5,778*u.K
R0 = const.R_sun

#Boltzmann and Wien's constants:
sg = const.sigma_sb
b = const.b_wien

#Scaling/ exponent values:
a = [3.5, 4]
eta = [0.57,0.8]
rf = 3.8/2.9

'''
    ------------------------------------------------------------------------------------------------------------
    This part is to store the directory information in a variable.
    The first part of this is to check if blender is returning a directory on it's own.
    The second part is to use the path of this script.
'''
dir = os.path.dirname(bpy.data.filepath)
if not dir or dir != "":
    dir = os.path.dirname(__file__)
if dir not in sys.path:
    sys.path.append(dir)
print(os.path.dirname(__file__))

    


def clean_slate():
    '''
    Simple function to clean up the scene,
    it deletes all objects and material from the scene
    '''
    for o in bpy.context.scene.objects:
        if o.type in ['MESH', 'EMPTY']:
            o.select_set(True)
        else:
            o.select_set(False)
    #
    bpy.ops.object.delete()
    bpy.data.lights.remove(bpy.data.lights[0])
    bpy.data.cameras.remove(bpy.data.cameras[0])
    #
    for m in bpy.data.materials:
        #add making env strength 0
        bpy.data.materials.remove(m)


class Star():

    def __init__(self, Mass, Radius, Time):
        self.M = Mass*u.M_sun
        self.rad = Radius
        self.t = Time
        if self.M > 1*u.M_sun:
            self.idx = 0
        else:
            self.idx = 1

    def props(self):

        if self.M<2*u.Msun:
            L = (L0*pow((self.M/(1*u.Msun)),a[1])).to(u.Lsun)
        else:
            L = (1.4*L0*pow((self.M/(1*u.Msun)),a[0])).to(u.L_sun)
            T = 0.9*T0*pow((self.M/(1*u.Msun)),0.625)
            R = rf*(R0*pow((self.M/(1*u.Msun)), eta[self.idx])).to(u.Rsun)
            peak_wavelength = (b/T).to(u.nm)
            lifetime = 10 *pow(self.M,-2.5)

        return [L, T, R, peak_wavelength, lifetime], [L.value, T.value, R.value, peak_wavelength.value, lifetime]

#The first list of values returned contains units and is an astropy Quantity object, while the second list contains just values
#Although radius has been calculated in the props function, it may not work well for all stars



    def coefflin(self):
        T=self.props[1,1]
        c1=-2.77873204e-09
        c2=-6.40472261e-05
        a1= pow(2.7, c1*T*T+c2*T)
        return a1
        
    def makeLimbDarkening(self):
        bpy.ops.mesh.primitive_uv_sphere_add(radius=1, enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(0.9, 0.9, 0.9))
        limb = bpy.context.active_object
        limb_mat = bpy.data.materials.new(name = "limb_darkening")
        limb.data.materials.append(limb_mat)
        #accessing nodes in material tab
        limb_mat.use_nodes = True
        
        limb_nodes = limb_mat.node_tree.nodes
        
        #removing node
        principled_bsdf = limb_nodes['Principled BSDF']
        limb_nodes.remove(principled_bsdf)

        #accessing material output node
        material_output = limb_nodes.get("Material Output")



        #new nodes
        geo = limb_nodes.new(type="ShaderNodeNewGeometry")
        vectra1 = limb_nodes.new(type="ShaderNodeVectorTransform")
        vectra1.convert_to = 'CAMERA'

        vectra2 = limb_nodes.new(type="ShaderNodeVectorTransform")
        vectra2.convert_to = 'CAMERA'
        
        dot_prod = limb_nodes.new( type = "ShaderNodeVectorMath")
        dot_prod.operation = 'DOT_PRODUCT'

        len1 = limb_nodes.new( type = "ShaderNodeVectorMath")
        len1.operation = 'LENGTH'

        len2 = limb_nodes.new( type = "ShaderNodeVectorMath")
        len2.operation = 'LENGTH'

        mult1 = limb_nodes.new( type = 'ShaderNodeMath')
        mult1.operation = 'MULTIPLY'
        mult1.inputs[1].default_value = 0.2

        mult2 = limb_nodes.new( type = 'ShaderNodeMath')
        mult2.operation = 'MULTIPLY'
        mult2.inputs[1].default_value = 0.2

        mult3 = limb_nodes.new( type = 'ShaderNodeMath')
        mult3.operation = 'MULTIPLY'

        div1 = limb_nodes.new( type = 'ShaderNodeMath')
        div1.operation = 'DIVIDE'

        cam_data = limb_nodes.new(type="ShaderNodeCameraData")

        div2 = limb_nodes.new( type = 'ShaderNodeMath')
        div2.operation = 'DIVIDE'
        div2.inputs[0].default_value = 1.8

        arccos = limb_nodes.new( type = 'ShaderNodeMath')
        arccos.operation = 'ARCCOSINE'

        sine1 = limb_nodes.new( type = 'ShaderNodeMath')
        sine1.operation = 'SINE'

        sine2 = limb_nodes.new( type = 'ShaderNodeMath')
        sine2.operation = 'SINE'

        div3 = limb_nodes.new( type = 'ShaderNodeMath')
        div3.operation = 'DIVIDE'

        power = limb_nodes.new( type = 'ShaderNodeMath')
        power.operation = 'POWER'
        power.inputs[1].default_value = 2.0

        sub1 = limb_nodes.new( type = 'ShaderNodeMath')
        sub1.operation = 'SUBTRACT'
        sub1.inputs[0].default_value = 1.0

        sqrt = limb_nodes.new( type = 'ShaderNodeMath')
        sqrt.operation = 'SQRT'

        mult4 = limb_nodes.new( type = 'ShaderNodeMath')
        mult4.operation = 'MULTIPLY'
        mult4.inputs[0].default_value = self.coefflin()

        sub2 = limb_nodes.new( type = 'ShaderNodeMath')
        sub2.operation = 'SUBTRACT'
        sub2.inputs[0].default_value = 1.0

        mult5 = limb_nodes.new( type = 'ShaderNodeMath')
        mult5.operation = 'MULTIPLY'
        mult5.inputs[0].default_value = 9.0

        gamma = limb_nodes.new( type = 'ShaderNodeGamma')
        gamma.inputs[0].default_value = (0.332452, 0.119538, 0.00560548, 1)
    
        

        #new nodes 
        br_con = limb_nodes.new( type = 'ShaderNodeBrightContrast')
        br_con.inputs[1].default_value = -0.560
        br_con.inputs[2].default_value = 0.200
        fresnel = limb_nodes.new(type = 'ShaderNodeFresnel')
        fresnel.inputs[0].default_value = 10.450
        emission = limb_nodes.new(type = 'ShaderNodeEmission')
        principled_bsdf = limb_nodes.new(type = 'ShaderNodeBsdfPrincipled')
        mix = limb_nodes.new(type = 'ShaderNodeMixShader')
        
        #adding nodes to invert colour
        sep_hsv = limb_nodes.new( type = 'ShaderNodeSeparateHSV')
        add1 = limb_nodes.new( type = 'ShaderNodeMath')
        add1.operation = 'ADD'
        add1.inputs[1].default_value = 0.5
        com_hsv = limb_nodes.new( type = 'ShaderNodeCombineHSV')
        sep_rgb = limb_nodes.new( type = 'ShaderNodeSeparateRGB')
        map_r1 = limb_nodes.new( type = 'ShaderNodeMapRange')
        map_r1.inputs[1].default_value = 0
        map_r1.inputs[2].default_value = 1
        map_r1.inputs[3].default_value = 1
        map_r1.inputs[4].default_value = 0
        map_r2 = limb_nodes.new( type = 'ShaderNodeMapRange')
        map_r2.inputs[1].default_value = 0
        map_r2.inputs[2].default_value = 1
        map_r2.inputs[3].default_value = 1
        map_r2.inputs[4].default_value = 0
        map_r3 = limb_nodes.new( type = 'ShaderNodeMapRange')
        map_r3.inputs[1].default_value = 0
        map_r3.inputs[2].default_value = 1
        map_r3.inputs[3].default_value = 1
        map_r3.inputs[4].default_value = 0

        com_rgb = limb_nodes.new( type = 'ShaderNodeCombineRGB')
        
        
        #making links
        limb_links = limb_mat.node_tree.links
        
        link1 = limb_links.new(geo.outputs[1], vectra1.inputs[0])
        link2 = limb_links.new(geo.outputs[4], vectra2.inputs[0])
        link3 = limb_links.new(geo.outputs[1], len1.inputs[0])
        link4 = limb_links.new(geo.outputs[4], len2.inputs[0])
        
        link5 = limb_links.new(vectra1.outputs[0], dot_prod.inputs[0])
        link6 = limb_links.new(vectra2.outputs[0], dot_prod.inputs[1])
        
        link7 = limb_links.new(len1.outputs[0], mult1.inputs[0])
        link8 = limb_links.new(len2.outputs[0], mult2.inputs[0])
        
        link9 = limb_links.new(dot_prod.outputs[1], div1.inputs[0])
        link10 = limb_links.new(dot_prod.outputs[1], gamma.inputs[1])
        
        link10b = limb_links.new(len1.outputs[1], mult1.inputs[0])
        link10c = limb_links.new(len2.outputs[1], mult2.inputs[0])
        
        link11 = limb_links.new(mult1.outputs[0], mult3.inputs[0])
        link12 = limb_links.new(mult2.outputs[0], mult3.inputs[1])
        
        link13 = limb_links.new(mult3.outputs[0], div1.inputs[1])
        
        link14 = limb_links.new(div1.outputs[0], arccos.inputs[0])
        link15 = limb_links.new(arccos.outputs[0], sine1.inputs[0])
        link16 = limb_links.new(sine1.outputs[0], div3.inputs[0])
        
        link17 = limb_links.new(cam_data.outputs[2], div2.inputs[1])
        link18 = limb_links.new(div2.outputs[0], sine2.inputs[0])
        link19 = limb_links.new(sine2.outputs[0], div3.inputs[1])
        
        link20 = limb_links.new(div3.outputs[0], power.inputs[0])
        link21 = limb_links.new(power.outputs[0], sub1.inputs[1])
        link22 = limb_links.new(sub1.outputs[0], sqrt.inputs[0])
        link23 = limb_links.new(sqrt.outputs[0], mult4.inputs[1])
        link24 = limb_links.new(mult4.outputs[0], sub2.inputs[1])
        link25 = limb_links.new(sub2.outputs[0], mult5.inputs[1])
        
        
        
        #linking nodes to invert colour
        link26 = limb_links.new(gamma.outputs[0], sep_hsv.inputs[0])
        link27 = limb_links.new(sep_hsv.outputs[0], add1.inputs[0])
        link28 = limb_links.new(add1.outputs[0], com_hsv.inputs[0])
        link29 = limb_links.new(sep_hsv.outputs[1], com_hsv.inputs[1])
        link30 = limb_links.new(sep_hsv.outputs[2], com_hsv.inputs[2])
        link31 = limb_links.new(com_hsv.outputs[0], sep_rgb.inputs[0])
        
        link32 = limb_links.new(sep_rgb.outputs[0], map_r1.inputs[0])
        link33 = limb_links.new(sep_rgb.outputs[1], map_r2.inputs[0])
        link34 = limb_links.new(sep_rgb.outputs[2], map_r3.inputs[0])
        
        link35 = limb_links.new(map_r1.outputs[0], com_rgb.inputs[0])
        link36 = limb_links.new(map_r2.outputs[0], com_rgb.inputs[1])
        link37 = limb_links.new(map_r3.outputs[0], com_rgb.inputs[2])
        
        link38 = limb_links.new(com_rgb.outputs[0], br_con.inputs[0])
        
        
        #---
        
        link39 = limb_links.new(br_con.outputs[0], emission.inputs[0])
        link40 = limb_links.new(mult5.outputs[0], emission.inputs[1])
        
        link41 = limb_links.new(fresnel.outputs[0], mix.inputs[0])
        link42 = limb_links.new(emission.outputs[0], mix.inputs[1])
        link43 = limb_links.new(principled_bsdf.outputs[0], mix.inputs[2])
        
        link44 = limb_links.new(mix.outputs[0], material_output.inputs[0])

        return limb


    def makeFlares(self):
        

        PhStarScale=1 #Scale of the place holder star
        FlareScale=0.2 #Scale of indidvidual flares
        FlareAmount=0.5 #controls the number of flares that can be seen on surface
        flareDensity= 25*self.props[1,1]/self.M*SCALE_FAC_MASS*SCALE_FAC_TEMP # controls how dense/thick the flares look, inversely prop to mass
        StreakColor=[1,0.089493,0,1] #color of flares
        StreakEmmissionStrength= 500/self.M*SCALE_FAC_MASS #strength of flares, inversely prop to mass




        #object for streaks
        bpy.ops.mesh.primitive_plane_add(enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))
        streak=bpy.context.active_object
        streak.scale=[.07,3.16,1]
        bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)

        # start material for streaks
        mat = bpy.data.materials.new(name="Flarmat")
        streak.data.materials.append(mat)
        mat.use_nodes = True
        nodes = mat.node_tree.nodes
        links = mat.node_tree.links

        mat_out = nodes.get("Material Output")
        bsdf = nodes.get("Principled BSDF")
        bsdf.inputs[17].default_value=StreakColor
        tex_coord = nodes.new(type="ShaderNodeTexCoord")
        map = nodes.new(type="ShaderNodeMapping")

        noise_tex1 = nodes.new(type="ShaderNodeTexNoise")
        noise_tex1.inputs[2].default_value = 6.1
        col_ramp = nodes.new(type="ShaderNodeValToRGB")
        col_ramp.color_ramp.elements[1].color = [0.216637,0.216637,0.216637,1]

        noise_tex2 = nodes.new(type="ShaderNodeTexNoise")
        noise_tex2.inputs[2].default_value = 29.4

        ramp1=nodes.new(type="ShaderNodeValToRGB")
        ramp1.color_ramp.elements[0].position = 0.34
        ramp1.color_ramp.elements[1].position = 0.44
        ramp1.color_ramp.elements.new(0.5)
        ramp1.color_ramp.elements.new(0.56)
        ramp1.color_ramp.elements.new(0.66)
        ramp1.color_ramp.elements[0].color=(0,0,0,1)
        ramp1.color_ramp.elements[1].color=(0.066875,0.066875,0.066875,1)
        ramp1.color_ramp.elements[2].color=(1,1,1,1)
        ramp1.color_ramp.elements[3].color=(0.066875,0.066875,0.066875,1)
        ramp1.color_ramp.elements[4].color=(0,0,0,1)

        sepXYZ=nodes.new(type="ShaderNodeSeparateXYZ")

        ramp2=nodes.new(type="ShaderNodeValToRGB")
        ramp2.color_ramp.interpolation='B_SPLINE'
        ramp2.color_ramp.elements[0].position = 0.11
        ramp2.color_ramp.elements[1].position = 0.23
        ramp2.color_ramp.elements.new(0.37)
        ramp2.color_ramp.elements.new(0.63)
        ramp2.color_ramp.elements.new(0.77)
        ramp2.color_ramp.elements.new(0.89)
        ramp2.color_ramp.elements[0].color=(0,0,0,1)
        ramp2.color_ramp.elements[1].color=(1,1,1,1)
        ramp2.color_ramp.elements[2].color=(0.075153,0.075153,0.075153,1)
        ramp2.color_ramp.elements[3].color=(0.075153,0.075153,0.075153,1)
        ramp2.color_ramp.elements[4].color=(1,1,1,1)
        ramp2.color_ramp.elements[5].color=(0,0,0,1)

        mix1=nodes.new(type="ShaderNodeMixRGB")
        mix1.blend_type='MULTIPLY'
        mix1.inputs[0].default_value=1
        mix2=nodes.new(type="ShaderNodeMixRGB")
        mix2.blend_type='MULTIPLY'
        mix2.inputs[0].default_value=1

        math=nodes.new(type="ShaderNodeMath")
        math.operation='MULTIPLY'
        math.inputs[1].default_value=StreakEmmissionStrength

        links.new(tex_coord.outputs[2], noise_tex1.inputs[1])
        links.new(tex_coord.outputs[2], map.inputs[0])
        links.new(tex_coord.outputs[2], sepXYZ.inputs[0])

        links.new(noise_tex1.outputs[1], col_ramp.inputs[0])
        links.new(col_ramp.outputs[0], map.inputs[1])

        links.new(sepXYZ.outputs[1], ramp2.inputs[0])
        links.new(map.outputs[0], ramp1.inputs[0])
        links.new(map.outputs[0], noise_tex2.inputs[1])

        links.new(ramp1.outputs[0], mix1.inputs[1])
        links.new(noise_tex2.outputs[0], mix1.inputs[2])

        links.new(ramp2.outputs[0], mix2.inputs[1])
        links.new(mix1.outputs[0], mix2.inputs[2])

        links.new(mix2.outputs[0], math.inputs[0])
        links.new(mix2.outputs[0], bsdf.inputs[0])
        links.new(mix2.outputs[0], bsdf.inputs[19])
        links.new(math.outputs[0], bsdf.inputs[18])

        #end material for streaks

        #start modifiers on streak
        bpy.ops.curve.primitive_bezier_circle_add(radius=1, enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))
        curve=bpy.context.active_object

        sub=streak.modifiers.new('subdiv','SUBSURF')
        sub.subdivision_type='SIMPLE'
        sub.levels=6

        cur=streak.modifiers.new('curve','CURVE')
        cur.object=curve
        cur.deform_axis='POS_Y'

        bpy.ops.object.select_all(action='DESELECT')
        bpy.context.view_layer.objects.active = streak
        bpy.ops.object.modifier_apply(modifier="subdiv")
        bpy.ops.object.modifier_apply(modifier="curve")

        dec=streak.modifiers.new('dec','DECIMATE')
        dec.decimate_type='DISSOLVE'
        dec.angle_limit=0.1*3.14/180
        bpy.ops.object.modifier_apply(modifier="dec")
        #end modifiers on streaks

        #flare object
        bpy.ops.mesh.primitive_circle_add(radius=0.115, fill_type='NGON', enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))
        flare=bpy.context.active_object

        #start geonodes for flares
        geoflare=flare.modifiers.new('geoflare','NODES')
        gfnodes=flare.modifiers["geoflare"].node_group.nodes
        gflinks=flare.modifiers["geoflare"].node_group.links

        inp=gfnodes.get("Group Input")
        outp=gfnodes.get("Group Output")
        dist=gfnodes.new(type="GeometryNodePointDistribute")
        dist.inputs[2].default_value=flareDensity

        randscale=gfnodes.new(type="GeometryNodeAttributeRandomize")
        randscale.inputs[1].default_value='scale'
        randscale.inputs[4].default_value=0.5
        randscale.inputs[5].default_value=1

        randrot=gfnodes.new(type="GeometryNodeAttributeRandomize")
        randrot.data_type='FLOAT_VECTOR'
        randrot.operation='ADD'
        randrot.inputs[1].default_value='rotation'
        randrot.inputs[2].default_value=[-0.7,0,0]
        randrot.inputs[3].default_value=[0.7,0,0]

        ptrot=gfnodes.new(type="GeometryNodeRotatePoints")
        ptrot.inputs[6].default_value=[-90*3.14/180,197.2*3.14/180,0]

        ptinst=gfnodes.new(type="GeometryNodePointInstance")
        ptinst.inputs[1].default_value=streak

        gflinks.new(inp.outputs[0],dist.inputs[0])
        gflinks.new(dist.outputs[0],randscale.inputs[0])
        gflinks.new(randscale.outputs[0],randrot.inputs[0])
        gflinks.new(randrot.outputs[0],ptrot.inputs[0])
        gflinks.new(ptrot.outputs[0],ptinst.inputs[0])
        gflinks.new(ptinst.outputs[0],outp.inputs[0])

        #end geo nodes for flares

        #star object(wont be visible) needs to be almost same size as star
        bpy.ops.mesh.primitive_uv_sphere_add(radius=1, enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(PhStarScale, PhStarScale, PhStarScale))
        star=bpy.context.active_object

        #start geonodes for star
        geostar=star.modifiers.new('geostar','NODES')
        gsnodes=star.modifiers["geostar"].node_group.nodes
        gslinks=star.modifiers["geostar"].node_group.links

        inp=gsnodes.get("Group Input")
        outp=gsnodes.get("Group Output")
        dist=gsnodes.new(type="GeometryNodePointDistribute")
        dist.inputs[2].default_value=FlareAmount

        randscale=gsnodes.new(type="GeometryNodeAttributeRandomize")
        randscale.inputs[1].default_value='scale'
        randscale.inputs[4].default_value=0
        randscale.inputs[5].default_value=FlareScale

        randrot=gsnodes.new(type="GeometryNodeAttributeRandomize")
        randrot.data_type='FLOAT_VECTOR'
        randrot.operation='ADD'
        randrot.inputs[1].default_value='rotation'
        randrot.inputs[2].default_value=[0,0,0]
        randrot.inputs[3].default_value=[0,0.4,0]


        ptinst=gsnodes.new(type="GeometryNodePointInstance")
        ptinst.inputs[1].default_value=flare

        gslinks.new(inp.outputs[0],dist.inputs[0])
        gslinks.new(dist.outputs[0],randscale.inputs[0])
        gslinks.new(randscale.outputs[0],randrot.inputs[0])
        gslinks.new(randrot.outputs[0],ptinst.inputs[0])
        gslinks.new(ptinst.outputs[0],outp.inputs[0])

        #end geonodes for star

        #hiding unnecessary stuff in Bin collection
        master=bpy.context.scene.view_layers['View Layer'].active_layer_collection.collection
        bin = bpy.data.collections.new("Bin")
        master.children.link(bin)

        master.objects.unlink(streak) 
        bin.objects.link(streak)
        master.objects.unlink(flare) 
        bin.objects.link(flare)
        master.objects.unlink(curve) 
        bin.objects.link(curve)

        bin.hide_render=True
        bin.hide_viewport=True

        return streak, flare, star, curve

    def makeMainSequenceStar(self):
        
        #bpy.context.scene.render.engine = 'BLENDER_EEVEE'
        #bpy.context.scene.eevee.use_bloom = True
        #bpy.context.scene.eevee.use_motion_blur = True

        bpy.ops.mesh.primitive_uv_sphere_add()
        ta=bpy.data.objects["Sphere"]
        bpy.ops.object.shade_smooth()
        ta.scale=(1,1,1) #scale of star
        ta.name="Main Sequence Star"
                
        ta.modifiers.new("part", 'PARTICLE_SYSTEM')
        part=ta.particle_systems[0]
        settings=part.settings

        sub=ta.modifiers.new("Subsurf", 'SUBSURF')

        #material for star
        tam=bpy.data.materials.new(name="base")
        ta.data.materials.append(tam)

        tam.use_nodes=True
        Nodes=tam.node_tree.nodes
        Links=tam.node_tree.links

        moutput=Nodes.get('Material Output')
        texc=Nodes.new(type='ShaderNodeTexCoord')
        map=Nodes.new(type="ShaderNodeMapping")

        mus=Nodes.new(type='ShaderNodeTexMusgrave')
        mus.musgrave_dimensions='4D'
        mus.musgrave_type='MULTIFRACTAL'
        mus.inputs[1].default_value=1.123
        mus.inputs[2].default_value=9.90
        mus.inputs[3].default_value=8.0
        mus.inputs[4].default_value=0.0
        mus.inputs[5].default_value=3.0

        noise=Nodes.new(type='ShaderNodeTexNoise')
        noise.inputs[2].default_value=25.4
        noise.inputs[3].default_value=13.5
        noise.inputs[4].default_value=0.5

        RGB=Nodes.new(type='ShaderNodeMixRGB')
        RGB.inputs[0].default_value=0.542
        cramp=Nodes.new(type='ShaderNodeValToRGB')
        cramp.color_ramp.interpolation='B_SPLINE'
        cramp.color_ramp.elements[0].position=0.42
        cramp.color_ramp.elements[1].position=0.505
        cramp.color_ramp.elements[0].color=(1.000,0.049,0.000,1.000) #might wanna consider temperature
        cramp.color_ramp.elements[1].color=(0.002,0.000,0.000117,1.000)

        ems=Nodes.new(type='ShaderNodeEmission')
        ems.inputs[1].default_value=self.props[1,1]/SCALE_FAC_TEMP

        Link1=Links.new(texc.outputs[2],map.inputs[0])
        Link2=Links.new(map.outputs[0],mus.inputs[0])
        Link3=Links.new(map.outputs[0],noise.inputs[0])
        Link4=Links.new(mus.outputs[0],RGB.inputs[1])
        Link5=Links.new(noise.outputs[1],RGB.inputs[2])
        Link6=Links.new(RGB.outputs[0],cramp.inputs[0])
        Link7=Links.new(cramp.outputs[0],ems.inputs[0])

        mix2=Nodes.new(type='ShaderNodeMixShader')
        mix2.inputs[0].default_value=0.542

        Link8=Links.new(ems.outputs[0],mix2.inputs[1])
        Link21=Links.new(mix2.outputs[0],moutput.inputs[0])

        #particle
        bpy.ops.mesh.primitive_uv_sphere_add()
        pa=bpy.data.objects["Sphere"]
        bpy.ops.object.shade_smooth()
        pa.name="Particle"
        pa.scale=(0.01,0.01,0.01)

        pam=bpy.data.materials.new(name="emitter")
        pa.data.materials.append(pam)

        pam.use_nodes=True
        Nodes1=pam.node_tree.nodes
        Links1=pam.node_tree.links

        moutput1=Nodes1.get('Material Output')
        ems1=Nodes1.new(type='ShaderNodeEmission')
        ems1.inputs[1].default_value=50
        ems1.inputs[0].default_value=(1.000,0.358,0.000,1.000)

        Link10=Links1.new(ems1.outputs[0],moutput1.inputs[0])

        #Emitter
        settings.type='EMITTER'
        settings.count=50000
        settings.frame_start=1
        settings.frame_end=50  ##
        settings.emit_from='FACE'
        settings.brownian_factor=0.13
        settings.drag_factor=0.042
        settings.render_type='OBJECT'
        settings.instance_object=bpy.data.objects["Particle"]

        bpy.ops.object.effector_add(type='GUIDE')
        cg=bpy.data.objects["CurveGuide"]
        cg.field.type='GUIDE'
        cg.field.guide_kink_type='CURL'
        cg.field.guide_kink_axis='Z'
        cg.field.guide_kink_frequency=1.5
        cg.field.guide_kink_amplitude=10
        cg.field.guide_kink_shape=0.9

        bpy.ops.object.effector_add(type='TURBULENCE')
        tb=bpy.data.objects["Turbulence"]
        tb.location=(1.01,0,0)
        tb.field.strength=10
        tb.field.size=7.38
        tb.field.seed=11

        #lava 
        texc1=Nodes.new(type='ShaderNodeTexCoord')
        map1=Nodes.new(type="ShaderNodeMapping")

        noise1=Nodes.new(type='ShaderNodeTexNoise')
        noise1.inputs[2].default_value=7.1
        noise1.inputs[3].default_value=16.0
        noise1.inputs[4].default_value=0.5
        noise1.inputs[5].default_value=1

        noise2=Nodes.new(type='ShaderNodeTexNoise')
        noise2.inputs[2].default_value=4.6
        noise2.inputs[3].default_value=16.0
        noise2.inputs[4].default_value=0.5
        noise2.inputs[5].default_value=1

        cramp1=Nodes.new(type='ShaderNodeValToRGB')
        cramp1.color_ramp.interpolation='LINEAR'
        cramp1.color_ramp.elements[0].position=0.5
        cramp1.color_ramp.elements[1].position=1.0

        prin=Nodes.new(type='ShaderNodeBsdfPrincipled')
        bump=Nodes.new(type='ShaderNodeBump')
        ems2=Nodes.new(type='ShaderNodeEmission')
        ems2.inputs[1].default_value=self.props[1,1]/(SCALE_FAC_TEMP*100)
        ems2.inputs[0].default_value=(1,0.665,0.263,1)

        mix=Nodes.new(type='ShaderNodeMixShader')

        Link11=Links.new(texc1.outputs[0],map1.inputs[0])
        Link12=Links.new(map1.outputs[0],noise1.inputs[0])
        Link13=Links.new(noise1.outputs[1],noise2.inputs[0])
        Link14=Links.new(noise2.outputs[0],cramp1.inputs[0])
        Link15=Links.new(cramp1.outputs[0],mix.inputs[0])
        Link16=Links.new(cramp1.outputs[0],bump.inputs[2])
        Link17=Links.new(bump.outputs[0],prin.inputs[20])
        Link18=Links.new(prin.outputs[0],mix.inputs[1])
        Link19=Links.new(ems2.outputs[0],mix.inputs[2])
        Link20=Links.new(mix.outputs[0],mix2.inputs[2])

        #world stars
        scene = bpy.context.scene
        space = bpy.data.worlds.new("space")
        scene.world = space
        space.use_nodes = True
        space_nodes = space.node_tree.nodes

        #accessing nodes
        space = space_nodes.get("Background")
        world1 = space_nodes.get("World Output")

        texc2=space_nodes.new(type='ShaderNodeTexCoord')
        map2=space_nodes.new(type="ShaderNodeMapping")

        noise3=space_nodes.new(type='ShaderNodeTexNoise')
        noise3.inputs[2].default_value=43.4
        noise3.inputs[3].default_value=3.6
        noise3.inputs[4].default_value=0.3
        noise3.inputs[5].default_value=0.3
        noise4=space_nodes.new(type='ShaderNodeTexNoise')
        noise4.inputs[3].default_value=2.0
        noise4.inputs[2].default_value=600

        cramp2=space_nodes.new(type='ShaderNodeValToRGB')
        cramp2.color_ramp.elements[0].position=0.695
        cramp3=space_nodes.new(type='ShaderNodeValToRGB')
        cramp3.color_ramp.elements[0].position=0.468
        math1=space_nodes.new(type='ShaderNodeMath')
        math1.inputs[1].default_value=500

        #space_links = space.node_tree.links

        #space_links.new(texc2.outputs[0],map2.inputs[0])
        #space_links.new(map2.outputs[0],noise3.inputs[0])
        #space_links.new(noise3.outputs[0],cramp3.inputs[0])
        #space_links.new(cramp3.outsputs[0],math1.inputs[0])
        #space_links.new(math1.outputs[0],back.inputs[1])
        #space_links.new(back.outputs[0],world1.inputs[0])
        #space_links.new(noise4.outputs[0],cramp2.inputs[0])
        #space_links.new(cramp2.outputs[0],back.inputs[0])

        return ta, pa, cg, tb
  

    def makeBlueGiant(self):

        bpy.ops.mesh.primitive_uv_sphere_add(radius=1, enter_editmode=False, align='WORLD', location=(0, 0, 0))
        #scaling to 10fm
        bpy.ops.transform.resize(value=(10, 10, 10), orient_type='GLOBAL', orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)), orient_matrix_type='GLOBAL', mirror=True, use_proportional_edit=False, proportional_edit_falloff='SMOOTH', proportional_size=1, use_proportional_connected=False, use_proportional_projected=False)

        #referencing sphere as earth
        bluegiant = bpy.context.active_object

        #shade smooth sphere
        mesh = bluegiant.data
        for face in mesh.polygons:
            face.use_smooth = True
            
            
        #adding base material
        bluegiant_mat = bpy.data.materials.new(name = "bluegiant_base")
        bluegiant.data.materials.append(bluegiant_mat)

        #accessing nodes in material tab
        bluegiant_mat.use_nodes = True
        bluegiant_nodes = bluegiant_mat.node_tree.nodes

        #removing node
        principled_bsdf = bluegiant_nodes['Principled BSDF']
        bluegiant_nodes.remove(principled_bsdf)

        #accessing material output node 
        material_output = bluegiant_nodes.get("Material Output")

        #adding new nodes
        emission = bluegiant_nodes.new(type = 'ShaderNodeEmission')
        tex_coord = bluegiant_nodes.new(type = 'ShaderNodeTexCoord')
        mapping = bluegiant_nodes.new(type = 'ShaderNodeMapping')
        colorramp1 = bluegiant_nodes.new(type = 'ShaderNodeValToRGB')
        colorramp2 = bluegiant_nodes.new(type = 'ShaderNodeValToRGB')
        colorramp3 = bluegiant_nodes.new(type = 'ShaderNodeValToRGB')
        noise1 = bluegiant_nodes.new(type = 'ShaderNodeTexNoise')
        noise2 = bluegiant_nodes.new(type = 'ShaderNodeTexNoise')
        noise3 = bluegiant_nodes.new(type = 'ShaderNodeTexNoise')
        multiply = bluegiant_nodes.new(type = 'ShaderNodeMixRGB')
        add = bluegiant_nodes.new(type = 'ShaderNodeMixRGB')
        color = bluegiant_nodes.new(type = 'ShaderNodeMixRGB')

        #changing node parameters
        mapping.vector_type = 'POINT'

        noise1.noise_dimensions = '4D'
        noise1.inputs[1].default_value = 0.0
        noise1.inputs[2].default_value = 30.0
        noise1.inputs[3].default_value = 15.8
        noise1.inputs[4].default_value = 0.825
        noise1.inputs[5].default_value = 0.3

        noise2.noise_dimensions = '4D'
        noise2.inputs[1].default_value = 0.0
        noise2.inputs[2].default_value = 30.0
        noise2.inputs[3].default_value = 15.8
        noise2.inputs[4].default_value = 0.825
        noise2.inputs[5].default_value = 0.3

        noise3.noise_dimensions = '4D'
        noise3.inputs[1].default_value = 0.2
        noise3.inputs[2].default_value = 8.40
        noise3.inputs[3].default_value = 1.50
        noise3.inputs[4].default_value = 0.783
        noise3.inputs[5].default_value = 0.660

        colorramp1.color_ramp.interpolation = 'B_SPLINE'
        colorramp1.color_ramp.elements[0].position = 0.532
        colorramp1.color_ramp.elements[1].position = 0.729
        colorramp1.color_ramp.elements.new(0.897)
        colorramp1.color_ramp.elements[0].color = (0.0,0.0,0.0,1.0)
        colorramp1.color_ramp.elements[1].color = (0.0,0.0,0.0,1.0)
        colorramp1.color_ramp.elements[2].color = (1.0,1.0,1.0,1.0)

        colorramp2.color_ramp.interpolation = 'B_SPLINE'
        colorramp2.color_ramp.elements[0].position = 0.319
        colorramp2.color_ramp.elements[1].position = 0.666
        colorramp2.color_ramp.elements.new(0.916)
        colorramp2.color_ramp.elements[0].color = (0.0,0.0,0.0,1.0)
        colorramp2.color_ramp.elements[1].color = (0.0,0.0,0.0,1.0)
        colorramp2.color_ramp.elements[2].color = (1.0,1.0,1.0,1.0)

        colorramp3.color_ramp.interpolation = 'B_SPLINE'
        colorramp3.color_ramp.elements[0].position = 0.322
        colorramp3.color_ramp.elements[1].position = 0.570
        colorramp3.color_ramp.elements.new(0.957)
        colorramp3.color_ramp.elements[0].color = (0.009,0.009,0.009,1.0)
        colorramp3.color_ramp.elements[1].color = (0.044,0.047,0.047,1.0)
        colorramp3.color_ramp.elements[2].color = (1.0,1.0,1.0,1.0)

        multiply.blend_type = 'MULTIPLY'
        multiply.inputs[0].default_value = 1.0

        add.blend_type = 'ADD'
        add.inputs[0].default_value = 1.0

        color.blend_type = 'COLOR'
        color.inputs[0].default_value = 1.0
        color.inputs[2].default_value[0] = 0.011
        color.inputs[2].default_value[1] = 0.129
        color.inputs[2].default_value[2] = 0.500

        if self.M < 10*M0:
            emission.inputs[1].default_value = 9*self.props[1,1]/SCALE_FAC_TEMP
        else:
            emission.inputs[1].default_value = 300 *self.props[1,0]

        #linking nodes
        bluegiant_links = bluegiant_mat.node_tree.links

        link1 = bluegiant_links.new(tex_coord.outputs[3], mapping.inputs[0])
        link2 = bluegiant_links.new(mapping.outputs[0], noise1.inputs[0])
        link3 = bluegiant_links.new(mapping.outputs[0], noise2.inputs[0])
        link4 = bluegiant_links.new(mapping.outputs[0], noise3.inputs[0])
        link5 = bluegiant_links.new(noise1.outputs[0], colorramp1.inputs[0])
        link6 = bluegiant_links.new(noise2.outputs[0], colorramp2.inputs[0])
        link7 = bluegiant_links.new(noise3.outputs[0], colorramp3.inputs[0])
        link8 = bluegiant_links.new(colorramp2.outputs[0], multiply.inputs[1])
        link9 = bluegiant_links.new(colorramp3.outputs[0], multiply.inputs[2])
        link10 = bluegiant_links.new(multiply.outputs[0], add.inputs[2])
        link11 = bluegiant_links.new(colorramp1.outputs[0], add.inputs[1])
        link12 = bluegiant_links.new(add.outputs[0], color.inputs[1])
        link13 = bluegiant_links.new(color.outputs[0], emission.inputs[0])
        link14 = bluegiant_links.new(emission.outputs[0], material_output.inputs[0])

        #compositing 
        bpy.data.scenes["Scene"].use_nodes = True
        comp= bpy.data.scenes["Scene"].node_tree
        comp_nodes = comp.nodes

        render_layer= comp_nodes.new(type = 'CompositorNodeRLayers')
        composite = comp_nodes.new(type = 'CompositorNodeComposite')
        glare = comp_nodes.new(type = 'CompositorNodeGlare')

        glare.glare_type = 'FOG_GLOW'
        glare.quality = 'HIGH'
        glare.mix = 0.0
        glare.threshold = 0.7
        glare.size = 9.0

        comp_links = comp.links
        link15 = comp_links.new(render_layer.outputs[0], glare.inputs[0])
        link16 = comp_links.new(glare.outputs[0], composite.inputs[0])

        return bluegiant


    def makeProtoStar(self):
        bpy.ops.mesh.primitive_uv_sphere_add(radius=1, enter_editmode=False, align='WORLD', location=(0, 0, 0))

        #scaling to 10fm
        bpy.ops.transform.resize(value=(0.5, 0.5, 0.5), orient_type='GLOBAL', orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)), orient_matrix_type='GLOBAL', mirror=True, use_proportional_edit=False, proportional_edit_falloff='SMOOTH', proportional_size=1, use_proportional_connected=False, use_proportional_projected=False)


        protostar = bpy.context.active_object
        #shade smooth sphere
        mesh = protostar.data
        for face in mesh.polygons:
            face.use_smooth = True
            

        #adding base material
        protostar_mat = bpy.data.materials.new(name = "protostar_base")
        protostar.data.materials.append(protostar_mat)

        #accessing nodes in material tab
        protostar_mat.use_nodes = True
        protostar_nodes = protostar_mat.node_tree.nodes

        #removing node
        principled_bsdf = protostar_nodes['Principled BSDF']
        protostar_nodes.remove(principled_bsdf)

        #accessing material output node 
        material_output = protostar_nodes.get("Material Output")

        #adding new nodes
        emission = protostar_nodes.new(type = 'ShaderNodeEmission')
        tex_coord = protostar_nodes.new(type = 'ShaderNodeTexCoord')
        mapping = protostar_nodes.new(type = 'ShaderNodeMapping')
        colorramp1 = protostar_nodes.new(type = 'ShaderNodeValToRGB')
        noise1 = protostar_nodes.new(type = 'ShaderNodeTexNoise')
        musgrave1 = protostar_nodes.new(type = 'ShaderNodeTexMusgrave')
        mix1 = protostar_nodes.new(type = 'ShaderNodeMixRGB')

        #changing node parameters
        mapping.vector_type = 'POINT'
        mapping.inputs[2].default_value[0] = 13.107 

        #texture coordinate
        tex_coord.object = protostar

        musgrave1.musgrave_dimensions ='4D'
        musgrave1.inputs[1].default_value = 1.123
        musgrave1.inputs[2].default_value = 9.0
        musgrave1.inputs[3].default_value = 8.0
        musgrave1.inputs[4].default_value = 0.0
        musgrave1.inputs[5].default_value = 3.0

        noise1.noise_dimensions = '3D'
        noise1.inputs[1].default_value = 25.0
        noise1.inputs[2].default_value = 11.0
        noise1.inputs[3].default_value = 0.5
        noise1.inputs[4].default_value = 0.0

        colorramp1.color_ramp.interpolation = 'B_SPLINE'
        colorramp1.color_ramp.elements[0].position = 0.417
        colorramp1.color_ramp.elements[1].position = 0.475
        colorramp1.color_ramp.elements[0].color = (0.019,0.153,0.416,1.0)
        colorramp1.color_ramp.elements[1].color = (0.0,0.0,0.0,1.0)

        emission.inputs[1].default_value = self.props[1,1]/SCALE_FAC_TEMP


        #linking nodes
        protostar_links = protostar_mat.node_tree.links

        link1 = protostar_links.new(tex_coord.outputs[2], mapping.inputs[0])
        link2 = protostar_links.new(mapping.outputs[0], musgrave1.inputs[0])
        link3 = protostar_links.new(mapping.outputs[0], noise1.inputs[0])
        link4 = protostar_links.new(musgrave1.outputs[0], mix1.inputs[1])
        link5 = protostar_links.new(noise1.outputs[1], mix1.inputs[2])
        link6 = protostar_links.new(mix1.outputs[0], colorramp1.inputs[0])
        link7 = protostar_links.new(colorramp1.outputs[0], emission.inputs[0])
        link8 = protostar_links.new(emission.outputs[0], material_output.inputs[0])




        #bpy.context.scene.eevee.use_bloom = True
        #bpy.context.scene.eevee.use_gtao = True


        #Adding Icosphere
        protostar.modifiers.new("part", 'PARTICLE_SYSTEM')
        part=protostar.particle_systems[0]
        settings=part.settings

        #particle

        bpy.ops.mesh.primitive_ico_sphere_add(enter_editmode=False, align='WORLD', location=(0, 60
        , 0), scale=(0.8, 0.8, 0.8 ))

        particle=bpy.data.objects["Icosphere"]
        particle.name="Particle"

        #particle =bpy.data.materials.new(name="emitter")
        #particle.data.materials.append(particle)


        #adding base material
        particle_mat = bpy.data.materials.new(name = "particle_base")
        particle.data.materials.append(particle_mat)

        #accessing nodes in material tab
        particle_mat.use_nodes = True
        particle_nodes = particle_mat.node_tree.nodes

        #removing node
        principled_bsdf = particle_nodes['Principled BSDF']
        particle_nodes.remove(principled_bsdf)

        #accessing material output node 
        material_output = particle_nodes.get("Material Output")

        #Adding nodes
        mixshader1 = particle_nodes.new(type = 'ShaderNodeMixShader')
        Emission2 = particle_nodes.new(type = 'ShaderNodeEmission')
        Transparent1 = particle_nodes.new(type = 'ShaderNodeBsdfTransparent')
        Colorramp2 = particle_nodes.new(type ='ShaderNodeValToRGB')
        Multiply2 = particle_nodes.new(type = 'ShaderNodeMath')
        Divide2 = particle_nodes.new(type = 'ShaderNodeMath')

        Particleinfo2 = particle_nodes.new(type = 'ShaderNodeParticleInfo')


        Multiply2.operation = 'MULTIPLY'

        Divide2.operation ='DIVIDE'

        Multiply2.inputs[1].default_value = 2.2
        Emission2.inputs[1].default_value = 10.0

        #Colorramp
        Colorramp2.color_ramp.interpolation = 'CONSTANT'
        Colorramp2.color_ramp.elements[0].position = 0.167
        Colorramp2.color_ramp.elements[1].position = 0.424
        Colorramp2.color_ramp.elements.new(0.906)
        Colorramp2.color_ramp.elements[0].color = (0.103,0.067,0.244,1.0)
        Colorramp2.color_ramp.elements[1].color = (0.733,0.084,0.088,1.0)
        Colorramp2.color_ramp.elements[2].color = (1.0,0.116,0.246,1.0)

        #creating links
        #linking nodes
        particle_links = particle_mat.node_tree.links

        link9 = particle_links.new(Particleinfo2.outputs[2], Divide2.inputs[0])
        link10 = particle_links.new(Particleinfo2.outputs[3], Divide2.inputs[1])
        link11 = particle_links.new(Divide2.outputs[0], Multiply2.inputs[0])
        link12 = particle_links.new(Divide2.outputs[0], Colorramp2.inputs[0])
        link13 = particle_links.new(Colorramp2.outputs[0], Emission2.inputs[0])
        link14 = particle_links.new(Multiply2.outputs[0], mixshader1.inputs[0])
        link15 = particle_links.new(Emission2.outputs[0], mixshader1.inputs[1])
        link16 = particle_links.new(Transparent1.outputs[0], mixshader1.inputs[2])
        link17 = particle_links.new(mixshader1.outputs[0], material_output.inputs[0])


        #Emitter
        settings.type='EMITTER'
        settings.count=10000
        settings.frame_end=100  
        settings.lifetime = 100
        settings.render_type='OBJECT'
        settings.effector_weights.gravity = 0
        settings.instance_object=bpy.data.objects["Particle"]


        #Adding turbulence
        bpy.ops.object.effector_add(type='TURBULENCE', enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(100, 100, 100))
        turbulence = bpy.data.objects["Turbulence"]
        bpy.context.object.field.strength = 50
        bpy.context.scene.render.frame_map_old = 50
        bpy.context.object.field.flow = 4
        bpy.context.object.field.noise = 10

        return protostar, particle, turbulence
    

    def makeWhiteDwarf(self):
        bpy.ops.mesh.primitive_uv_sphere_add(radius=1, enter_editmode=False, align='WORLD', location=(0, 0, 0))

        #scaling to 10fm
        bpy.ops.transform.resize(value=(0.5, 0.5, 0.5), orient_type='GLOBAL', orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)), orient_matrix_type='GLOBAL', mirror=True, use_proportional_edit=False, proportional_edit_falloff='SMOOTH', proportional_size=1, use_proportional_connected=False, use_proportional_projected=False)


        protostar = bpy.context.active_object
        #shade smooth sphere
        mesh = protostar.data
        for face in mesh.polygons:
            face.use_smooth = True
            

        #adding base material
        protostar_mat = bpy.data.materials.new(name = "protostar_base")
        protostar.data.materials.append(protostar_mat)

        #accessing nodes in material tab
        protostar_mat.use_nodes = True
        protostar_nodes = protostar_mat.node_tree.nodes

        #removing node
        principled_bsdf = protostar_nodes['Principled BSDF']
        protostar_nodes.remove(principled_bsdf)

        #accessing material output node 
        material_output = protostar_nodes.get("Material Output")

        #adding new nodes
        emission = protostar_nodes.new(type = 'ShaderNodeEmission')
        tex_coord = protostar_nodes.new(type = 'ShaderNodeTexCoord')
        mapping = protostar_nodes.new(type = 'ShaderNodeMapping')
        colorramp1 = protostar_nodes.new(type = 'ShaderNodeValToRGB')
        noise1 = protostar_nodes.new(type = 'ShaderNodeTexNoise')
        musgrave1 = protostar_nodes.new(type = 'ShaderNodeTexMusgrave')
        mix1 = protostar_nodes.new(type = 'ShaderNodeMixRGB')

        #changing node parameters
        mapping.vector_type = 'POINT'
        mapping.inputs[2].default_value[0] = 13.107 

        #texture coordinate
        tex_coord.object = protostar

        musgrave1.musgrave_dimensions ='4D'
        musgrave1.inputs[1].default_value = 1.123
        musgrave1.inputs[2].default_value = 9.0
        musgrave1.inputs[3].default_value = 8.0
        musgrave1.inputs[4].default_value = 0.0
        musgrave1.inputs[5].default_value = 3.0

        noise1.noise_dimensions = '3D'
        noise1.inputs[1].default_value = 25.0
        noise1.inputs[2].default_value = 11.0
        noise1.inputs[3].default_value = 0.5
        noise1.inputs[4].default_value = 0.0

        colorramp1.color_ramp.interpolation = 'B_SPLINE'
        colorramp1.color_ramp.elements[0].position = 0.417
        colorramp1.color_ramp.elements[1].position = 0.475
        colorramp1.color_ramp.elements[0].color = (0.019,0.153,0.416,1.0)
        colorramp1.color_ramp.elements[1].color = (0.0,0.0,0.0,1.0)

        emission.inputs[1].default_value = self.props[1,0]/1000


        #linking nodes
        protostar_links = protostar_mat.node_tree.links

        link1 = protostar_links.new(tex_coord.outputs[2], mapping.inputs[0])
        link2 = protostar_links.new(mapping.outputs[0], musgrave1.inputs[0])
        link3 = protostar_links.new(mapping.outputs[0], noise1.inputs[0])
        link4 = protostar_links.new(musgrave1.outputs[0], mix1.inputs[1])
        link5 = protostar_links.new(noise1.outputs[1], mix1.inputs[2])
        link6 = protostar_links.new(mix1.outputs[0], colorramp1.inputs[0])
        link7 = protostar_links.new(colorramp1.outputs[0], emission.inputs[0])
        link8 = protostar_links.new(emission.outputs[0], material_output.inputs[0])




        #bpy.context.scene.eevee.use_bloom = True
        #bpy.context.scene.eevee.use_gtao = True






        #bpy.context.scene.frame_start = 250
        #bpy.context.scene.frame_end = 10

        return protostar
    

    def makeRedGiant(self):
        bpy.ops.mesh.primitive_uv_sphere_add(radius=1, enter_editmode=False, align='WORLD', location=(0, 0, 0))
        #scaling to 10fm
        bpy.ops.transform.resize(value=(10, 10, 10), orient_type='GLOBAL', orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)), orient_matrix_type='GLOBAL', mirror=True, use_proportional_edit=False, proportional_edit_falloff='SMOOTH', proportional_size=1, use_proportional_connected=False, use_proportional_projected=False)

        #referencing sphere as earth
        redgiant = bpy.context.active_object

        #shade smooth sphere
        mesh = redgiant.data
        for face in mesh.polygons:
            face.use_smooth = True
            
            
        #adding base material
        redgiant_mat = bpy.data.materials.new(name = "redgiant_base")
        redgiant.data.materials.append(redgiant_mat)

        #accessing nodes in material tab
        redgiant_mat.use_nodes = True
        redgiant_nodes = redgiant_mat.node_tree.nodes

        #removing node
        r_principled_bsdf = redgiant_nodes['Principled BSDF']
        redgiant_nodes.remove(r_principled_bsdf)

        #accessing material output node 
        r_material_output = redgiant_nodes.get("Material Output")

        #adding new nodes
        r_emission = redgiant_nodes.new(type = 'ShaderNodeEmission')
        r_tex_coord = redgiant_nodes.new(type = 'ShaderNodeTexCoord')
        r_mapping = redgiant_nodes.new(type = 'ShaderNodeMapping')
        r_colorramp1 = redgiant_nodes.new(type = 'ShaderNodeValToRGB')
        r_colorramp2 = redgiant_nodes.new(type = 'ShaderNodeValToRGB')
        r_colorramp3 = redgiant_nodes.new(type = 'ShaderNodeValToRGB')
        r_noise1 = redgiant_nodes.new(type = 'ShaderNodeTexNoise')
        r_noise2 = redgiant_nodes.new(type = 'ShaderNodeTexNoise')
        r_noise3 = redgiant_nodes.new(type = 'ShaderNodeTexNoise')
        r_multiply = redgiant_nodes.new(type = 'ShaderNodeMixRGB')
        r_add = redgiant_nodes.new(type = 'ShaderNodeMixRGB')
        r_color1 = redgiant_nodes.new(type = 'ShaderNodeMixRGB')
        r_color2 = redgiant_nodes.new(type = 'ShaderNodeMixRGB')

        #changing node parameters
        r_mapping.vector_type = 'POINT'

        r_noise1.noise_dimensions = '4D'
        r_noise1.inputs[1].default_value = 0.0
        r_noise1.inputs[2].default_value = 25.0
        r_noise1.inputs[3].default_value = 15.8
        r_noise1.inputs[4].default_value = 0.825
        r_noise1.inputs[5].default_value = 0.3

        r_noise2.noise_dimensions = '4D'
        r_noise2.inputs[1].default_value = 0.0
        r_noise2.inputs[2].default_value = 20.0
        r_noise2.inputs[3].default_value = 15.8
        r_noise2.inputs[4].default_value = 0.825
        r_noise2.inputs[5].default_value = 0.3

        r_noise3.noise_dimensions = '4D'
        r_noise3.inputs[1].default_value = 0.2
        r_noise3.inputs[2].default_value = 5.0
        r_noise3.inputs[3].default_value = 1.50
        r_noise3.inputs[4].default_value = 0.783
        r_noise3.inputs[5].default_value = 0.660

        r_colorramp1.color_ramp.interpolation = 'B_SPLINE'
        r_colorramp1.color_ramp.elements[0].position = 0.486
        r_colorramp1.color_ramp.elements[1].position = 0.615
        r_colorramp1.color_ramp.elements.new(0.944)
        r_colorramp1.color_ramp.elements[0].color = (0.0,0.0,0.0,1.0)
        r_colorramp1.color_ramp.elements[1].color = (0.0,0.0,0.0,1.0)
        r_colorramp1.color_ramp.elements[2].color = (1.0,1.0,1.0,1.0)

        r_colorramp2.color_ramp.interpolation = 'B_SPLINE'
        r_colorramp2.color_ramp.elements[0].position = 0.275
        r_colorramp2.color_ramp.elements[1].position = 0.705
        r_colorramp2.color_ramp.elements.new(0.822)
        r_colorramp2.color_ramp.elements[0].color = (0.0,0.0,0.0,1.0)
        r_colorramp2.color_ramp.elements[1].color = (0.0,0.0,0.0,1.0)
        r_colorramp2.color_ramp.elements[2].color = (1.0,1.0,1.0,1.0)

        r_colorramp3.color_ramp.interpolation = 'B_SPLINE'
        r_colorramp3.color_ramp.elements[0].position = 0.387
        r_colorramp3.color_ramp.elements[1].position = 0.770
        r_colorramp3.color_ramp.elements.new(0.917)
        r_colorramp3.color_ramp.elements[0].color = (0.009,0.009,0.009,1.0)
        r_colorramp3.color_ramp.elements[1].color = (0.044,0.047,0.047,1.0)
        r_colorramp3.color_ramp.elements[2].color = (1.0,1.0,1.0,1.0)

        r_multiply.blend_type = 'MULTIPLY'
        r_multiply.inputs[0].default_value = 1.0

        r_add.blend_type = 'ADD'
        r_add.inputs[0].default_value = 1.0

        r_color1.blend_type = 'COLOR'
        r_color1.inputs[0].default_value = 1.0
        r_color1.inputs[2].default_value[0] = 0.5
        r_color1.inputs[2].default_value[1] = 0.006179
        r_color1.inputs[2].default_value[2] = 0.001444

        r_color2.blend_type = 'COLOR'
        r_color2.inputs[0].default_value = 0.135
        r_color2.inputs[2].default_value[0] = 0.5
        r_color2.inputs[2].default_value[1] = 0.072221
        r_color2.inputs[2].default_value[2] = 0.003184

        if self.M < 10*M0:
            r_emission.inputs[1].default_value = 9*self.props[1,1]/SCALE_FAC_TEMP
        else:
            r_emission.inputs[1].default_value = 300 *self.props[1,0]

        #linking nodes
        redgiant_links = redgiant_mat.node_tree.links

        link1 = redgiant_links.new(r_tex_coord.outputs[3], r_mapping.inputs[0])
        link2 = redgiant_links.new(r_mapping.outputs[0], r_noise1.inputs[0])
        link3 = redgiant_links.new(r_mapping.outputs[0], r_noise2.inputs[0])
        link4 = redgiant_links.new(r_mapping.outputs[0], r_noise3.inputs[0])
        link5 = redgiant_links.new(r_noise1.outputs[0], r_colorramp1.inputs[0])
        link6 = redgiant_links.new(r_noise2.outputs[0], r_colorramp2.inputs[0])
        link7 = redgiant_links.new(r_noise3.outputs[0], r_colorramp3.inputs[0])
        link8 = redgiant_links.new(r_colorramp2.outputs[0], r_multiply.inputs[1])
        link9 = redgiant_links.new(r_colorramp3.outputs[0], r_multiply.inputs[2])
        link10 = redgiant_links.new(r_multiply.outputs[0], r_add.inputs[2])
        link11 = redgiant_links.new(r_colorramp1.outputs[0], r_add.inputs[1])
        link12 = redgiant_links.new(r_add.outputs[0], r_color1.inputs[1])
        link13 = redgiant_links.new(r_color1.outputs[0], r_color2.inputs[1])
        link14 = redgiant_links.new(r_color2.outputs[0], r_emission.inputs[0])
        link15 = redgiant_links.new(r_emission.outputs[0], r_material_output.inputs[0])

        #compositing 
        bpy.data.scenes["Scene"].use_nodes = True
        r_comp= bpy.data.scenes["Scene"].node_tree
        r_comp_nodes = r_comp.nodes

        r_render_layer= r_comp_nodes.get("Render Layers")
        r_composite = r_comp_nodes.get("Composite")
        r_glare = r_comp_nodes.new(type = 'CompositorNodeGlare')

        r_glare.glare_type = 'FOG_GLOW'
        r_glare.quality = 'LOW'
        r_glare.mix = 0.2
        r_glare.threshold = 1.0
        r_glare.size = 8.0

        r_comp_links = r_comp.links
        link16 = r_comp_links.new(r_render_layer.outputs[0], r_glare.inputs[0])
        link17 = r_comp_links.new(r_glare.outputs[0], r_composite.inputs[0])

        return redgiant

    def makeCamera():
        camera_data = bpy.data.cameras.new(name='Camera')
        camera = bpy.data.objects.new('Camera', camera_data)
        camera.location = (19.1811, -23.4414, 7.12684)
        bpy.context.scene.camera = bpy.data.objects["Camera"]
        print(camera.location)
        return camera


if __name__ == "__main__":

    argv = sys.argv

    if "--" not in argv:
        argv = []
    else:
        argv = argv[argv.index("--") + 1:]

        usage_text = (
            "Run Blender to get an astrophysically accurate Star:"
            "  blender --background --python " + __file__ + " -- [options]"
        )

        parser = argparse.ArgumentParser(description=usage_text)
        parser.add_argument("--render-anim", "-ra", action='store_true', default=False, help="Use this flag if you want to render the default animation in the background.", dest="renderAnim")
        args = parser.parse_args(argv)

        
        
        if args.renderAnim:
            inputFromUser = [float(item) for item in input("Enter the parameters for the required star in the following order: 1. Mass (in solar mass units) 2. Radius (in solar radius units) 3. Time(in billion years) ")]
            star = Star(inputFromUser.split(','))
            #cam = star
            #cam.makeCamera()



            if star.t > star.props[1,4] and star.M < 5*const.M_sun:
                #make camera
                cam_rg = star
                cam_rg.makeCamera()
                #make red giant
                a = star
                a.makeRedGiant()
                c = star
                d = c.coefflin()
                c.makeLimbDarkening()
                b = star
                b.makeFlares()
            
            elif star.t > star.props[1,4] and star.M > 5*const.M_sun:
                #make camera
                cam_bg = star
                cam_bg.makeCamera()
                #blue giant
                a = star
                a.makeBlueGiant()
                c = star
                d = c.coefflin()
                c.makeLimbDarkening()
                b = star
                b.makeFlares()
            
        
            elif star.t == 0:
                #make camera
                cam_ps = star
                cam_ps.makeCamera()
                #proto star
                a = star
                a.makeProtoStar()
                c = star
                d = c.coefflin()
                c.makeLimbDarkening()
                b = star
                b.makeFlares()

            elif star.t > 1.1*star.props[1,4]:
                #make camera
                cam_wd = star
                cam_wd.makeCamera()
                #white dwarf
                a = star
                a.makeWhiteDwarf()
                c = star
                d = c.coefflin()
                c.makeLimbDarkening()
                b = star
                b.makeFlares()


            elif star.t < star.props[1,4]:
                #make camera
                cam_ms = star
                cam_ms.makeCamera()
                #main sequence star
                a = star
                a.makeMainSequenceStar()
                c = star
                d = c.coefflin()
                c.makeLimbDarkening()
                b = star
                b.makeFlares()
            
            bpy.ops.wm.save_as_mainfile(filepath="./starrender.blend")

            start_frame = 1
            bpy.context.scene.frame_start = start_frame
            end_frame = 20
            bpy.context.scene.frame_end = end_frame
            bpy.context.scene.render.fps = 24
            bpy.context.scene.render.engine = 'CYCLES'

            #bpy.context.scene.render.filepath = os.path.join(dir, 'Star')

            bpy.context.scene.render.image_settings.file_format = "FFMPEG"
            bpy.context.scene.cycles.samples = 16

            bpy.context.scene.render.filepath = os.path.join(__file__, 'StarRender')
            bpy.ops.render.render(animation=True, write_still=1)
