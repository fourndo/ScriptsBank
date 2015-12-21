# Info: Current History file: C:\Users\dominiquef\Dropbox\DOM_Projects\Meshing_scripts\Synthetic_tiles.prj\history-this-session-dominiquef-01e3428a-b195c731-f2cc-f7b3a8796340.gsc 


# ** ---------------------------------------------------------------- ** 
# ** New Session 
# ** Start Date: Thu Oct 31 17:15:20 2013

# ** User ID: dominiquef
gocad archive_new overwrite_existing true path C:/Users/dominiquef/Dropbox/DOM_Projects/Meshing_scripts\Synthetic_tiles.prj cs_name Default xy_unit_name m time_unit_name ms depth_unit_name m z_positive elevation x_positive east default_cs depth
# 
# 
gocad new_sgrid name "Tile1" origin 0. 0. 0. step_u 10 0. 0. step_v 0. 10 0. step_w 0. 0. -10 nu 10 nv 10 nw 10 cell_centered "true" coordinate_system_name "Default_depth" 
# 
# 
gocad on Style SGrid.Tile1  set attribute "*volume" value "true"
# 
# 
gocad on Style SGrid.Tile1 update
# 
# 
gocad update all cameras
# 
# 
gocad on GObj Tile1 create_property  property_name prop no_datavalue_specified 1 no_datavalue -99999 number_of_elements 1 property_kind "unknown" storage Memory
# 
# 
gocad on GObj "Tile1" set_value_property property "prop" value 1 region "everywhere" 
# 
# 
gocad on Style "Tile1" set_attribute attribute "*painted" value "on" 
# 
# 
gocad on Style "Tile1" set_attribute attribute "*painted*variable" value "prop" 
# 
# 
gocad on Style "Tile1" update 
# 
# 
gocad on Camera "DefaultCamera" set_camera_colormap object_name "Tile1" property_name "prop" show_colormap "false" 
# 
# 
gocad update_cameras_showing gobjs "Tile1" 
# 
# 
gocad copy gobj "Tile1" name "Tile2" copy_properties "true" copy_regions "false" copy_style "true" copy_points "true" copy_constraints "false" 
# 
# 
gocad copy gobj "Tile1" name "Tile3" copy_properties "true" copy_regions "false" copy_style "true" copy_points "true" copy_constraints "false" 
# 
# 
gocad gobjs_move objects "Tile2" translate "true" translation 20 20 0.0 rotate "false" origin 0.0 0.0 0.0 axis 0.0 0.0 1.0 angle 0 
# 
# 
gocad gobjs_move objects "Tile2" translate "true" translation 20 20 0.0 rotate "false" origin 0.0 0.0 0.0 axis 0.0 0.0 1.0 angle 0 
# 
# 
gocad gobjs_move objects "Tile3" translate "true" translation 20 -20 -10 rotate "false" origin 0.0 0.0 0.0 axis 0.0 0.0 1.0 angle 0 
# 
# 
gocad gobjs_move objects "Tile3" translate "true" translation 20 -20 -10 rotate "false" origin 0.0 0.0 0.0 axis 0.0 0.0 1.0 angle 0 
# 
# 
gocad gobjs_move objects "Tile3" translate "true" translation 0 10 0 rotate "false" origin 0.0 0.0 0.0 axis 0.0 0.0 1.0 angle 0 
# 
# 
gocad gobjs_move objects "Tile3" translate "true" translation 0 10 0 rotate "false" origin 0.0 0.0 0.0 axis 0.0 0.0 1.0 angle 0 
# 
# 
gocad on GObj "Tile2" set_value_property property "prop" value 2 region "everywhere" 
# 
# 
gocad on GObj "Tile3" set_value_property property "prop" value 3 region "everywhere" 
# 
# 
gocad on Style Tile2/properties/prop prop update
# 
# 
gocad update all cameras
# 
# 
gocad on Camera "Camera#0" set_top_view 
# 
# 
gocad on Camera "Camera#0" set_top_view 
# 
# 
gocad mirapotentialcreatemeshfilefromgridobj void_param "" name_ "Tile1" Mesh_File_name "C:\Users\dominiquef\Dropbox\DOM_Projects\Meshing_scripts\Test\Tile1.msh" 
# 
# 
gocad mirapotentialcreatemeshfilefromgridobj void_param "" name_ "Tile2" Mesh_File_name "C:\Users\dominiquef\Dropbox\DOM_Projects\Meshing_scripts\Test\Tile2.msh" 
# 
# 
gocad mirapotentialcreatemeshfilefromgridobj void_param "" name_ "Tile3" Mesh_File_name "C:\Users\dominiquef\Dropbox\DOM_Projects\Meshing_scripts\Test\Tile3.msh" 
# 
# 
gocad mirapotentialcreatemodelfilefromgridobj void_param "" name_ "Tile1" model_or_lower_bound_property "prop" bounds_file "false" upper_bound_property "prop" Model_File_name "C:\Users\dominiquef\Dropbox\DOM_Projects\Meshing_scripts\Test\Tile1.den" 
# 
# 
gocad mirapotentialcreatemodelfilefromgridobj void_param "" name_ "Tile2" model_or_lower_bound_property "prop" bounds_file "false" upper_bound_property "prop" Model_File_name "C:\Users\dominiquef\Dropbox\DOM_Projects\Meshing_scripts\Test\Tile2.den" 
# 
# 
gocad mirapotentialcreatemodelfilefromgridobj void_param "" name_ "Tile3" model_or_lower_bound_property "prop" bounds_file "false" upper_bound_property "prop" Model_File_name "C:\Users\dominiquef\Dropbox\DOM_Projects\Meshing_scripts\Test\Tile3.den" 
# Info: Saving style_catalog... 
# Info: Saving Tile1... 
# Info: Saving Tile2... 
# Info: Saving Tile3... 
# Info: Saving application_objects... 
# Info: Saving geobaselib... 
# Info: Saving geologic_features... 
