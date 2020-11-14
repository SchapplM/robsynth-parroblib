% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RRPRR8V2G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3*3x13]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 21:04
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR8V2G1A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:03:21
% EndTime: 2020-08-06 21:03:24
% DurationCPUTime: 2.87s
% Computational Cost: add. (4012->284), mult. (5220->564), div. (840->18), fcn. (4434->35), ass. (0->248)
t847 = 2 * pkin(1);
t709 = cos(qJ(2,3));
t678 = t709 * pkin(2);
t660 = t678 + pkin(1);
t711 = cos(qJ(2,2));
t679 = t711 * pkin(2);
t661 = t679 + pkin(1);
t713 = cos(qJ(2,1));
t680 = t713 * pkin(2);
t662 = t680 + pkin(1);
t659 = cos(pkin(7)) * pkin(3) + pkin(2);
t703 = sin(qJ(2,3));
t827 = sin(pkin(7)) * pkin(3);
t626 = t703 * t659 + t709 * t827;
t629 = t659 * t709 - t703 * t827;
t846 = t626 ^ 2 / t629 ^ 2;
t705 = sin(qJ(2,2));
t627 = t705 * t659 + t711 * t827;
t630 = t659 * t711 - t705 * t827;
t845 = t627 ^ 2 / t630 ^ 2;
t707 = sin(qJ(2,1));
t628 = t707 * t659 + t713 * t827;
t631 = t659 * t713 - t707 * t827;
t844 = t628 ^ 2 / t631 ^ 2;
t702 = pkin(5) + qJ(3,1);
t683 = -pkin(6) - t702;
t676 = 0.1e1 / t683;
t798 = t676 * t713;
t835 = -0.2e1 * pkin(2);
t755 = t798 * t835;
t831 = pkin(1) * t676;
t843 = -t755 / 0.2e1 + t831;
t701 = pkin(5) + qJ(3,2);
t682 = -pkin(6) - t701;
t674 = 0.1e1 / t682;
t803 = t674 * t711;
t756 = t803 * t835;
t832 = pkin(1) * t674;
t842 = -t756 / 0.2e1 + t832;
t700 = pkin(5) + qJ(3,3);
t681 = -pkin(6) - t700;
t672 = 0.1e1 / t681;
t808 = t672 * t709;
t757 = t808 * t835;
t833 = pkin(1) * t672;
t841 = -t757 / 0.2e1 + t833;
t792 = pkin(1) ^ 2 + pkin(5) ^ 2;
t834 = 2 * pkin(5);
t658 = (t834 + qJ(3,1)) * qJ(3,1) + t792;
t720 = pkin(2) ^ 2;
t765 = t676 * t713 ^ 2 * t720;
t840 = pkin(1) * t755 - t658 * t676 - t765;
t657 = (t834 + qJ(3,2)) * qJ(3,2) + t792;
t767 = t674 * t711 ^ 2 * t720;
t839 = pkin(1) * t756 - t657 * t674 - t767;
t656 = (t834 + qJ(3,3)) * qJ(3,3) + t792;
t769 = t672 * t709 ^ 2 * t720;
t838 = pkin(1) * t757 - t656 * t672 - t769;
t837 = -2 * pkin(1);
t687 = qJ(2,3) + pkin(7);
t830 = pkin(3) * cos(t687);
t688 = qJ(2,2) + pkin(7);
t829 = pkin(3) * cos(t688);
t689 = qJ(2,1) + pkin(7);
t828 = pkin(3) * cos(t689);
t653 = t678 + t830;
t650 = 0.1e1 / t653;
t715 = 0.2e1 * qJ(2,3);
t719 = pkin(3) ^ 2;
t793 = 0.2e1 * pkin(2) * pkin(3);
t826 = (sin(t715 + pkin(7)) * t793 + t720 * sin(t715) + t719 * sin(0.2e1 * t687) + (sin(t687) * pkin(3) + pkin(2) * t703) * t847) * t650;
t654 = t679 + t829;
t651 = 0.1e1 / t654;
t716 = 0.2e1 * qJ(2,2);
t825 = (sin(t716 + pkin(7)) * t793 + t720 * sin(t716) + t719 * sin(0.2e1 * t688) + (sin(t688) * pkin(3) + pkin(2) * t705) * t847) * t651;
t655 = t680 + t828;
t652 = 0.1e1 / t655;
t717 = 0.2e1 * qJ(2,1);
t824 = (sin(t717 + pkin(7)) * t793 + t720 * sin(t717) + t719 * sin(0.2e1 * t689) + (sin(t689) * pkin(3) + pkin(2) * t707) * t847) * t652;
t823 = t626 / t629;
t822 = t627 / t630;
t821 = t628 / t631;
t697 = legFrame(3,3);
t666 = sin(t697);
t669 = cos(t697);
t704 = sin(qJ(1,3));
t710 = cos(qJ(1,3));
t644 = -t666 * t704 + t669 * t710;
t820 = t644 * t672;
t698 = legFrame(2,3);
t667 = sin(t698);
t670 = cos(t698);
t706 = sin(qJ(1,2));
t712 = cos(qJ(1,2));
t645 = -t667 * t706 + t670 * t712;
t819 = t645 * t674;
t699 = legFrame(1,3);
t668 = sin(t699);
t671 = cos(t699);
t708 = sin(qJ(1,1));
t714 = cos(qJ(1,1));
t646 = -t668 * t708 + t671 * t714;
t818 = t646 * t676;
t647 = t666 * t710 + t669 * t704;
t817 = t647 * t672;
t648 = t667 * t712 + t670 * t706;
t816 = t648 * t674;
t649 = t668 * t714 + t671 * t708;
t815 = t649 * t676;
t814 = t650 * t703;
t813 = t651 * t705;
t812 = t652 * t707;
t673 = 0.1e1 / t681 ^ 2;
t690 = t703 ^ 2;
t807 = t673 * t690;
t806 = t673 * t700;
t805 = t673 * t703;
t804 = t673 * t709;
t675 = 0.1e1 / t682 ^ 2;
t691 = t705 ^ 2;
t802 = t675 * t691;
t801 = t675 * t701;
t800 = t675 * t705;
t799 = t675 * t711;
t677 = 0.1e1 / t683 ^ 2;
t692 = t707 ^ 2;
t797 = t677 * t692;
t796 = t677 * t702;
t795 = t677 * t707;
t794 = t677 * t713;
t790 = t673 * t846;
t789 = t675 * t845;
t788 = t677 * t844;
t787 = t673 * t823;
t786 = t675 * t822;
t785 = t677 * t821;
t784 = t672 * t823;
t783 = t674 * t822;
t782 = t676 * t821;
t781 = t644 * t673 * t647;
t780 = t645 * t675 * t648;
t779 = t646 * t677 * t649;
t778 = t672 * t814;
t777 = t650 * t808;
t776 = t700 * t814;
t775 = t674 * t813;
t774 = t651 * t803;
t773 = t701 * t813;
t772 = t676 * t812;
t771 = t652 * t798;
t770 = t702 * t812;
t768 = t703 * t804;
t766 = t705 * t799;
t764 = t707 * t794;
t760 = t826 / 0.2e1;
t759 = t825 / 0.2e1;
t758 = t824 / 0.2e1;
t754 = t703 * t784;
t753 = t709 * t784;
t752 = t705 * t783;
t751 = t711 * t783;
t750 = t707 * t782;
t749 = t713 * t782;
t748 = t709 * t781;
t747 = t644 * t787;
t746 = t711 * t780;
t745 = t645 * t786;
t744 = t713 * t779;
t743 = t646 * t785;
t742 = t647 * t787;
t741 = t648 * t786;
t740 = t649 * t785;
t739 = t672 * t776;
t738 = t674 * t773;
t737 = t676 * t770;
t736 = 0.2e1 * t700 * t787;
t735 = 0.2e1 * t701 * t786;
t734 = 0.2e1 * t702 * t785;
t733 = t768 * t823;
t732 = t766 * t822;
t731 = t764 * t821;
t730 = t650 * t754;
t729 = t651 * t752;
t728 = t652 * t750;
t727 = t644 * t778 + t645 * t775 + t646 * t772;
t726 = t644 * t777 + t645 * t774 + t646 * t771;
t725 = t647 * t778 + t648 * t775 + t649 * t772;
t724 = t647 * t777 + t648 * t774 + t649 * t771;
t723 = t728 + t729 + t730;
t722 = t650 * t753 + t651 * t751 + t652 * t749;
t643 = t649 ^ 2;
t642 = t648 ^ 2;
t641 = t647 ^ 2;
t640 = t646 ^ 2;
t639 = t645 ^ 2;
t638 = t644 ^ 2;
t637 = t661 * t712 - t706 * t682;
t636 = t660 * t710 - t704 * t681;
t635 = t714 * t662 - t708 * t683;
t634 = t708 * t662 + t714 * t683;
t633 = t706 * t661 + t712 * t682;
t632 = t704 * t660 + t710 * t681;
t613 = -pkin(2) * t812 - 0.2e1 * t702 * t782;
t612 = -pkin(2) * t813 - 0.2e1 * t701 * t783;
t611 = -pkin(2) * t814 - 0.2e1 * t700 * t784;
t610 = -t713 * pkin(5) * t652 + t750 * t847;
t609 = -t711 * pkin(5) * t651 + t752 * t847;
t608 = -t709 * pkin(5) * t650 + t754 * t847;
t607 = -pkin(5) * t812 + t749 * t837;
t606 = -pkin(5) * t813 + t751 * t837;
t605 = -pkin(5) * t814 + t753 * t837;
t604 = t634 * t671 + t635 * t668 + t649 * t828;
t603 = t633 * t670 + t637 * t667 + t648 * t829;
t602 = t632 * t669 + t636 * t666 + t647 * t830;
t601 = -t634 * t668 + t635 * t671 + t646 * t828;
t600 = -t633 * t667 + t637 * t670 + t645 * t829;
t599 = -t632 * t666 + t636 * t669 + t644 * t830;
t594 = t779 + t780 + t781;
t593 = (-t649 * t662 + t604) * t676;
t592 = (-t648 * t661 + t603) * t674;
t591 = (-t647 * t660 + t602) * t672;
t590 = (-t646 * t662 + t601) * t676;
t589 = (-t645 * t661 + t600) * t674;
t588 = (-t644 * t660 + t599) * t672;
t587 = t690 * t781 + t691 * t780 + t692 * t779;
t586 = 0.2e1 * t700 * t781 + 0.2e1 * t701 * t780 + 0.2e1 * t702 * t779;
t585 = (t744 + t746 + t748) * t847;
t584 = (-t703 * t781 - t705 * t780 - t707 * t779) * t847;
t583 = 0.2e1 * t703 * t748 + 0.2e1 * t705 * t746 + 0.2e1 * t707 * t744;
t582 = (-t662 * t821 + t758) * t676;
t581 = (-t661 * t822 + t759) * t674;
t580 = (-t660 * t823 + t760) * t672;
t579 = t740 + t741 + t742;
t578 = t743 + t745 + t747;
t577 = t690 * t742 + t691 * t741 + t692 * t740;
t576 = t690 * t747 + t691 * t745 + t692 * t743;
t575 = 0.2e1 * t647 * t733 + 0.2e1 * t648 * t732 + 0.2e1 * t649 * t731;
t574 = 0.2e1 * t644 * t733 + 0.2e1 * t645 * t732 + 0.2e1 * t646 * t731;
t573 = t843 * t604 + t840 * t649;
t572 = t843 * t601 + t840 * t646;
t571 = t842 * t603 + t839 * t648;
t570 = t842 * t600 + t839 * t645;
t569 = t841 * t602 + t838 * t647;
t568 = t841 * t599 + t838 * t644;
t567 = -t765 * t821 + (pkin(1) * t821 - t824 / 0.4e1) * t755 - pkin(2) * t770 - t658 * t782 + t758 * t831;
t566 = -t767 * t822 + (pkin(1) * t822 - t825 / 0.4e1) * t756 - pkin(2) * t773 - t657 * t783 + t759 * t832;
t565 = -t769 * t823 + (pkin(1) * t823 - t826 / 0.4e1) * t757 - pkin(2) * t776 - t656 * t784 + t760 * t833;
t1 = [t638 * t673 + t639 * t675 + t640 * t677, 0, 0, t638 * t807 + t639 * t802 + t640 * t797, 0.2e1 * t638 * t768 + 0.2e1 * t639 * t766 + 0.2e1 * t640 * t764, 0, 0, 0, (t638 * t804 + t639 * t799 + t640 * t794) * t847, (-t638 * t805 - t639 * t800 - t640 * t795) * t847, 0.2e1 * t638 * t806 + 0.2e1 * t639 * t801 + 0.2e1 * t640 * t796, -(t572 * t646 - t590 * t601) * t676 - (t570 * t645 - t589 * t600) * t674 - (t568 * t644 - t588 * t599) * t672, 1; t594, 0, 0, t587, t583, 0, 0, 0, t585, t584, t586, -(t573 * t646 - t593 * t601) * t676 - (t571 * t645 - t592 * t600) * t674 - (t569 * t644 - t591 * t599) * t672, 0; t578, 0, 0, t576, t574, -t727, -t726, 0, -t605 * t820 - t606 * t819 - t607 * t818, -t608 * t820 - t609 * t819 - t610 * t818, -t611 * t820 - t612 * t819 - t613 * t818, -(t567 * t646 - t582 * t601) * t676 - (t566 * t645 - t581 * t600) * t674 - (t565 * t644 - t580 * t599) * t672, 0; t594, 0, 0, t587, t583, 0, 0, 0, t585, t584, t586, -(t572 * t649 - t590 * t604) * t676 - (t570 * t648 - t589 * t603) * t674 - (t568 * t647 - t588 * t602) * t672, 0; t641 * t673 + t642 * t675 + t643 * t677, 0, 0, t641 * t807 + t642 * t802 + t643 * t797, 0.2e1 * t641 * t768 + 0.2e1 * t642 * t766 + 0.2e1 * t643 * t764, 0, 0, 0, (t641 * t804 + t642 * t799 + t643 * t794) * t847, (-t641 * t805 - t642 * t800 - t643 * t795) * t847, 0.2e1 * t641 * t806 + 0.2e1 * t642 * t801 + 0.2e1 * t643 * t796, -(t573 * t649 - t593 * t604) * t676 - (t571 * t648 - t592 * t603) * t674 - (t569 * t647 - t591 * t602) * t672, 1; t579, 0, 0, t577, t575, -t725, -t724, 0, -t605 * t817 - t606 * t816 - t607 * t815, -t608 * t817 - t609 * t816 - t610 * t815, -t611 * t817 - t612 * t816 - t613 * t815, -(t567 * t649 - t582 * t604) * t676 - (t566 * t648 - t581 * t603) * t674 - (t565 * t647 - t580 * t602) * t672, 0; t578, 0, 0, t576, t574, -t727, -t726, 0, t727 * pkin(5) + (t709 * t747 + t711 * t745 + t713 * t743) * t847, t726 * pkin(5) + (-t703 * t747 - t705 * t745 - t707 * t743) * t847, pkin(2) * t727 + t644 * t736 + t645 * t735 + t646 * t734, -(t572 * t821 - t590 * t758) * t676 - (t570 * t822 - t589 * t759) * t674 - (t568 * t823 - t588 * t760) * t672 + (t644 * t739 + t645 * t738 + t646 * t737) * pkin(2), 0; t579, 0, 0, t577, t575, -t725, -t724, 0, t725 * pkin(5) + (t709 * t742 + t711 * t741 + t713 * t740) * t847, t724 * pkin(5) + (-t703 * t742 - t705 * t741 - t707 * t740) * t847, pkin(2) * t725 + t647 * t736 + t648 * t735 + t649 * t734, -(t573 * t821 - t593 * t758) * t676 - (t571 * t822 - t592 * t759) * t674 - (t569 * t823 - t591 * t760) * t672 + (t647 * t739 + t648 * t738 + t649 * t737) * pkin(2), 0; t788 + t789 + t790, 0, 0, t690 * t790 + t691 * t789 + t692 * t788, 0.2e1 * t764 * t844 + 0.2e1 * t766 * t845 + 0.2e1 * t768 * t846, -0.2e1 * t723, -0.2e1 * t722, 0.1e1 / t655 ^ 2 + 0.1e1 / t654 ^ 2 + 0.1e1 / t653 ^ 2, pkin(5) * t723 - t605 * t784 - t606 * t783 - t607 * t782, pkin(5) * t722 - t608 * t784 - t609 * t783 - t610 * t782, pkin(2) * t723 - t611 * t784 - t612 * t783 - t613 * t782, -(t567 * t821 - t582 * t758) * t676 - (t566 * t822 - t581 * t759) * t674 - (t565 * t823 - t580 * t760) * t672 + (t700 * t730 + t701 * t729 + t702 * t728 + (t650 ^ 2 + t651 ^ 2 + t652 ^ 2) * pkin(2)) * pkin(2), 1;];
tau_reg  = t1;
