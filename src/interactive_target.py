#!/usr/bin/python

import rospy
import copy

from interactive_markers.interactive_marker_server import *
from interactive_markers.menu_handler import *
from visualization_msgs.msg import *
from geometry_msgs.msg import Point
from tf.broadcaster import TransformBroadcaster


# def frameCallback(msg):
#     global counter, br
#     time = rospy.Time.now()
#     br.sendTransform((0, 0, sin(counter / 140.0) * 2.0), (
#         0, 0, 0, 1.0), time, "base_link", "moving_frame")
#     counter += 1

def makeBox( msg ):
    marker = Marker()

    marker.type = Marker.SPHERE
    marker.scale.x = msg.scale * 0.02
    marker.scale.y = msg.scale * 0.02
    marker.scale.z = msg.scale * 0.02
    marker.color.r = 0.0
    marker.color.g = 1.0
    marker.color.b = 0.0
    marker.color.a = 1.0

    return marker

def makeBoxControl( msg ):
    control =  InteractiveMarkerControl()
    control.always_visible = True
    control.markers.append( makeBox(msg) )
    msg.controls.append( control )
    return control

def processFeedback(feedback, br):
    s = "Feedback from marker '" + feedback.marker_name
    s += "' / control '" + feedback.control_name + "'"

    mp = ""
    if feedback.mouse_point_valid:
        mp = " at " + str(feedback.mouse_point.x)
        mp += ", " + str(feedback.mouse_point.y)
        mp += ", " + str(feedback.mouse_point.z)
        mp += " in frame " + feedback.header.frame_id

    if feedback.event_type == InteractiveMarkerFeedback.BUTTON_CLICK:
        rospy.loginfo(s + ": button click" + mp + ".")
    elif feedback.event_type == InteractiveMarkerFeedback.MENU_SELECT:
        rospy.loginfo(s + ": menu item " + str(
            feedback.menu_entry_id) + " clicked" + mp + ".")
    elif feedback.event_type == InteractiveMarkerFeedback.POSE_UPDATE:
        rospy.loginfo(s + ": pose changed")
        br.sendTransform((feedback.pose.position.x,
                         feedback.pose.position.y,
                         feedback.pose.position.z),
                         (0, 0, 0, 1),
                         rospy.Time.now(),
                         "ds_target",
                         "world")
# TODO
#          << "\nposition = "
#          << feedback.pose.position.x
#          << ", " << feedback.pose.position.y
#          << ", " << feedback.pose.position.z
#          << "\norientation = "
#          << feedback.pose.orientation.w
#          << ", " << feedback.pose.orientation.x
#          << ", " << feedback.pose.orientation.y
#          << ", " << feedback.pose.orientation.z
#          << "\nframe: " << feedback.header.frame_id
#          << " time: " << feedback.header.stamp.sec << "sec, "
#          << feedback.header.stamp.nsec << " nsec" )
    elif feedback.event_type == InteractiveMarkerFeedback.MOUSE_DOWN:
        rospy.loginfo(s + ": mouse down" + mp + ".")
    elif feedback.event_type == InteractiveMarkerFeedback.MOUSE_UP:
        rospy.loginfo(s + ": mouse up" + mp + ".")
    server.applyChanges()
    

if __name__ == "__main__":
    rospy.init_node("target_interaction")
    br = TransformBroadcaster()

    # create a timer to update the published transforms
    # rospy.Timer(rospy.Duration(0.01), frameCallback)

    server = InteractiveMarkerServer("target_interaction")
    menu_handler = MenuHandler()
    pf_wrap = lambda fb: processFeedback(fb, br)
    
    menu_handler.insert("First Entry",
                        callback=pf_wrap)
    menu_handler.insert("Second Entry",
                        callback=pf_wrap)
    sub_menu_handle = menu_handler.insert("Submenu")
    menu_handler.insert(
        "First Entry", parent=sub_menu_handle, callback=pf_wrap)
    menu_handler.insert(
        "Second Entry", parent=sub_menu_handle, callback=pf_wrap)

    position = Point(0.5, 0.5, 0.5)

    int_marker = InteractiveMarker()
    int_marker.header.frame_id = "world"
    int_marker.pose.position = position
    int_marker.scale = 1

    int_marker.name = "motion target"
    int_marker.description = "Target point for DS motion generator."

    # insert a box
    makeBoxControl(int_marker)
    fixed = False
    # control = InteractiveMarkerControl()
    # control.orientation.w = 1
    # control.orientation.x = 1
    # control.orientation.y = 0
    # control.orientation.z = 0
    # control.name = "rotate_x"
    # control.interaction_mode = InteractiveMarkerControl.ROTATE_AXIS
    # if fixed:
    #     control.orientation_mode = InteractiveMarkerControl.FIXED
    # int_marker.controls.append(control)

    control = InteractiveMarkerControl()
    control.orientation.w = 1
    control.orientation.x = 1
    control.orientation.y = 0
    control.orientation.z = 0
    control.name = "move_x"
    control.interaction_mode = InteractiveMarkerControl.MOVE_AXIS
    if fixed:
        control.orientation_mode = InteractiveMarkerControl.FIXED
    int_marker.controls.append(control)

    # control = InteractiveMarkerControl()
    # control.orientation.w = 1
    # control.orientation.x = 0
    # control.orientation.y = 1
    # control.orientation.z = 0
    # control.name = "rotate_z"
    # control.interaction_mode = InteractiveMarkerControl.ROTATE_AXIS
    # if fixed:
    #     control.orientation_mode = InteractiveMarkerControl.FIXED
    # int_marker.controls.append(control)

    control = InteractiveMarkerControl()
    control.orientation.w = 1
    control.orientation.x = 0
    control.orientation.y = 1
    control.orientation.z = 0
    control.name = "move_z"
    control.interaction_mode = InteractiveMarkerControl.MOVE_AXIS
    if fixed:
        control.orientation_mode = InteractiveMarkerControl.FIXED
    int_marker.controls.append(control)

    # control = InteractiveMarkerControl()
    # control.orientation.w = 1
    # control.orientation.x = 0
    # control.orientation.y = 0
    # control.orientation.z = 1
    # control.name = "rotate_y"
    # control.interaction_mode = InteractiveMarkerControl.ROTATE_AXIS
    # if fixed:
    #     control.orientation_mode = InteractiveMarkerControl.FIXED
    # int_marker.controls.append(control)

    control = InteractiveMarkerControl()
    control.orientation.w = 1
    control.orientation.x = 0
    control.orientation.y = 0
    control.orientation.z = 1
    control.name = "move_y"
    control.interaction_mode = InteractiveMarkerControl.MOVE_AXIS
    if fixed:
        control.orientation_mode = InteractiveMarkerControl.FIXED
    int_marker.controls.append(control)

    server.insert(int_marker, pf_wrap)
    menu_handler.apply(server, int_marker.name)
    # make6DofMarker(
    #     False, InteractiveMarkerControl.MOVE_ROTATE_3D, position, True)

    server.applyChanges()

    rospy.spin()
