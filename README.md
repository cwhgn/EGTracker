An Equalised Global Graphical Model-Based Approach for Multi-Camera Object Tracking
===================

--------------------------

Weihua Chen, Lijun Cao, Xiaotang Chen, Kaiqi Huang 

-----------------------

National laboratory of pattern recognition, Chinese Academy of Sciences

--------------
Multi-camera non-overlapping visual object tracking system typically consists of two tasks: single camera object tracking and inter-camera object tracking. Since the state-of-the-art approaches are yet not perform perfectly in real scenes, the errors in single camera object tracking module would propagate into the module of inter-camera object tracking, resulting much lower overall performance. In order to address this problem, we develop an approach that jointly optimise the single camera object tracking and inter-camera object tracking in an equalised global graphical model. Such an approach has the advantage of guaranteeing a good overall tracking performance even when there are limited amount of false tracking in single camera object tracking. Besides, the similarity metrics used in our approach improve the compatibility of the metrics used in the two different tasks. Results show that our approach achieve the state-of-the-art results in multi-camera non-overlapping tracking datasets.


arXiv preprint:
http://arxiv.org/abs/1502.03532

demo:
http://youtu.be/GZ2u2tvzgi4
,
http://v.youku.com/v_show/id_XODkzMDgxNjM2.html

Code Prerequisites:
VS2008
Opencv1.0 + Opencv2.0 (The combination of opencv1.0 and 2.0 can be find in opencv file.)
FYI, as our project import some other libs, which use either opencv1.0 or opencv 2.0, so we have to require both opencv1.0 and opencv 2.0 in our project for these libs.



