% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RPRRR12V1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x14]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:21
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPRRR12V1G1A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:21:35
% EndTime: 2020-08-06 18:21:38
% DurationCPUTime: 2.86s
% Computational Cost: add. (10334->283), mult. (16140->557), div. (1761->18), fcn. (11982->18), ass. (0->244)
t633 = sin(qJ(3,3));
t592 = pkin(3) * t633 + qJ(2,3);
t582 = 0.1e1 / t592 ^ 2;
t639 = cos(qJ(3,3));
t645 = xDP(3);
t762 = t639 * t645;
t630 = legFrame(3,3);
t599 = sin(t630);
t602 = cos(t630);
t634 = sin(qJ(1,3));
t640 = cos(qJ(1,3));
t646 = xDP(2);
t647 = xDP(1);
t697 = t634 * t646 + t640 * t647;
t698 = -t634 * t647 + t640 * t646;
t560 = t698 * t599 + t697 * t602;
t614 = 0.1e1 / t633;
t785 = t560 * t614;
t648 = pkin(1) + pkin(5);
t613 = pkin(6) + t648;
t786 = t560 * t613;
t581 = 0.1e1 / t592;
t787 = t560 * t581;
t590 = t613 * t647;
t575 = -qJ(2,3) * t646 + t590;
t591 = t646 * t613;
t578 = qJ(2,3) * t647 + t591;
t539 = ((t575 * t640 + t578 * t634) * t602 + (-t575 * t634 + t578 * t640) * t599) * t633 + qJ(2,3) * t762 + (t633 * t762 + (-t697 * t599 + t698 * t602) * (t639 - 0.1e1) * (t639 + 0.1e1)) * pkin(3);
t798 = t539 * t614;
t583 = t581 * t582;
t799 = t539 * t583;
t524 = (0.2e1 * t582 * t762 - t799) * t785 - (-t786 + t798) * t582 * t787;
t664 = t633 ^ 2;
t615 = 0.1e1 / t664;
t629 = t645 ^ 2;
t766 = t629 / pkin(3) ^ 2;
t724 = t648 * t766;
t713 = t615 * t724;
t739 = t539 * t785;
t722 = 0.2e1 * t739;
t809 = 0.2e1 * qJ(2,3);
t817 = t524 * t809 + t582 * t722 + t713;
t635 = sin(qJ(3,2));
t593 = pkin(3) * t635 + qJ(2,2);
t585 = 0.1e1 / t593 ^ 2;
t641 = cos(qJ(3,2));
t761 = t641 * t645;
t631 = legFrame(2,3);
t600 = sin(t631);
t603 = cos(t631);
t636 = sin(qJ(1,2));
t642 = cos(qJ(1,2));
t695 = t636 * t646 + t642 * t647;
t696 = -t636 * t647 + t642 * t646;
t561 = t696 * t600 + t695 * t603;
t618 = 0.1e1 / t635;
t782 = t561 * t618;
t783 = t561 * t613;
t584 = 0.1e1 / t593;
t784 = t561 * t584;
t576 = -qJ(2,2) * t646 + t590;
t579 = qJ(2,2) * t647 + t591;
t540 = ((t576 * t642 + t579 * t636) * t603 + (-t576 * t636 + t579 * t642) * t600) * t635 + qJ(2,2) * t761 + (t635 * t761 + (-t695 * t600 + t696 * t603) * (t641 - 0.1e1) * (t641 + 0.1e1)) * pkin(3);
t796 = t540 * t618;
t586 = t584 * t585;
t797 = t540 * t586;
t525 = (0.2e1 * t585 * t761 - t797) * t782 - (-t783 + t796) * t585 * t784;
t667 = t635 ^ 2;
t619 = 0.1e1 / t667;
t711 = t619 * t724;
t738 = t540 * t782;
t721 = 0.2e1 * t738;
t810 = 0.2e1 * qJ(2,2);
t816 = t525 * t810 + t585 * t721 + t711;
t637 = sin(qJ(3,1));
t594 = -t637 * pkin(3) - qJ(2,1);
t588 = 0.1e1 / t594 ^ 2;
t643 = cos(qJ(3,1));
t760 = t645 * t643;
t632 = legFrame(1,3);
t601 = sin(t632);
t604 = cos(t632);
t638 = sin(qJ(1,1));
t644 = cos(qJ(1,1));
t693 = t638 * t646 + t644 * t647;
t694 = -t638 * t647 + t644 * t646;
t562 = t694 * t601 + t693 * t604;
t622 = 0.1e1 / t637;
t779 = t562 * t622;
t780 = t562 * t613;
t587 = 0.1e1 / t594;
t781 = t562 * t587;
t577 = -qJ(2,1) * t646 + t590;
t580 = qJ(2,1) * t647 + t591;
t541 = ((t577 * t644 + t580 * t638) * t604 + (-t577 * t638 + t580 * t644) * t601) * t637 + qJ(2,1) * t760 + (t637 * t760 + (-t693 * t601 + t694 * t604) * (t643 - 0.1e1) * (t643 + 0.1e1)) * pkin(3);
t794 = t541 * t622;
t589 = t587 * t588;
t795 = t541 * t589;
t526 = (0.2e1 * t588 * t760 + t795) * t779 + (-t780 + t794) * t588 * t781;
t670 = t637 ^ 2;
t623 = 0.1e1 / t670;
t709 = t623 * t724;
t737 = t541 * t779;
t720 = 0.2e1 * t737;
t811 = 0.2e1 * qJ(2,1);
t815 = t526 * t811 + t588 * t720 + t709;
t626 = t639 ^ 2;
t814 = t614 * t626;
t627 = t641 ^ 2;
t813 = t618 * t627;
t628 = t643 ^ 2;
t812 = t622 * t628;
t649 = qJ(2,3) ^ 2;
t652 = pkin(3) ^ 2;
t653 = 0.1e1 / pkin(3);
t759 = t645 * t653;
t692 = pkin(3) * t613 * t759;
t743 = qJ(2,3) * t759;
t752 = -t613 ^ 2 - t652;
t765 = t633 * t639;
t679 = -((-t614 * t692 * t765 + (t613 * t633 * t798 + ((t626 * t652 - t649 + t752) * t633 + (t626 - 0.1e1) * t809 * pkin(3)) * t560) * t581) * t582 + t613 * t799) * t785 + ((t581 * t639 * t786 + t614 * t645) * t633 + t614 * t743) * t581 * t615 * t645;
t650 = qJ(2,2) ^ 2;
t745 = qJ(2,2) * t759;
t764 = t635 * t641;
t678 = -((-t618 * t692 * t764 + (t613 * t635 * t796 + ((t627 * t652 - t650 + t752) * t635 + (t627 - 0.1e1) * t810 * pkin(3)) * t561) * t584) * t585 + t613 * t797) * t782 + ((t584 * t641 * t783 + t618 * t645) * t635 + t618 * t745) * t584 * t619 * t645;
t651 = qJ(2,1) ^ 2;
t747 = qJ(2,1) * t759;
t763 = t637 * t643;
t677 = (-(-t622 * t692 * t763 - (t613 * t637 * t794 + ((t628 * t652 - t651 + t752) * t637 + (t628 - 0.1e1) * t811 * pkin(3)) * t562) * t587) * t588 + t613 * t795) * t779 - ((-t587 * t643 * t780 + t622 * t645) * t637 + t622 * t747) * t587 * t623 * t645;
t808 = pkin(1) * t524;
t807 = pkin(1) * t525;
t806 = pkin(1) * t526;
t805 = qJ(2,1) * t589;
t804 = qJ(2,2) * t586;
t803 = qJ(2,3) * t583;
t802 = t524 * t581;
t801 = t525 * t584;
t800 = t526 * t587;
t557 = t560 ^ 2;
t793 = t557 * t582;
t792 = t557 * t583;
t558 = t561 ^ 2;
t791 = t558 * t585;
t790 = t558 * t586;
t559 = t562 ^ 2;
t789 = t559 * t588;
t788 = t559 * t589;
t778 = t582 * t639;
t777 = t585 * t641;
t776 = t588 * t643;
t775 = t614 * t639;
t774 = t615 * t814;
t773 = 0.1e1 / t664 ^ 2 * t639;
t772 = t618 * t641;
t771 = t619 * t813;
t770 = 0.1e1 / t667 ^ 2 * t641;
t769 = t622 * t643;
t768 = t623 * t812;
t767 = 0.1e1 / t670 ^ 2 * t643;
t702 = t743 * t787;
t758 = 0.2e1 * t702 + (-t713 + t817) * t639;
t712 = t766 * t774;
t757 = t817 * t633 + t648 * t712 - 0.2e1 * t702 * t775;
t703 = t745 * t784;
t756 = 0.2e1 * t703 + (-t711 + t816) * t641;
t710 = t766 * t771;
t755 = t816 * t635 + t648 * t710 - 0.2e1 * t703 * t772;
t704 = t747 * t781;
t754 = -0.2e1 * t704 + (-t709 + t815) * t643;
t708 = t766 * t768;
t753 = t815 * t637 + t648 * t708 + 0.2e1 * t704 * t769;
t748 = qJ(2,1) * t789;
t746 = qJ(2,2) * t791;
t744 = qJ(2,3) * t793;
t742 = t626 * t802;
t741 = t627 * t801;
t740 = t628 * t800;
t736 = t557 * t778;
t735 = t558 * t777;
t734 = t559 * t776;
t733 = t560 * t778;
t732 = t561 * t777;
t731 = t562 * t776;
t730 = t582 * (0.2e1 * t626 - 0.1e1) * t614;
t729 = t585 * (0.2e1 * t627 - 0.1e1) * t618;
t728 = t588 * (0.2e1 * t628 - 0.1e1) * t622;
t727 = t615 * t766;
t726 = t619 * t766;
t725 = t623 * t766;
t723 = 0.2e1 * t759;
t719 = t765 * t802;
t718 = t764 * t801;
t717 = t763 * t800;
t716 = t560 * t730;
t715 = t561 * t729;
t714 = t562 * t728;
t707 = t526 * t648 - t677 + t748;
t706 = t525 * t648 - t678 + t746;
t705 = t524 * t648 - t679 + t744;
t701 = (-t614 - t774) * t581;
t700 = (-t618 - t771) * t584;
t699 = (-t622 - t768) * t587;
t548 = -t727 - t793;
t691 = (t548 + t727) * t639;
t549 = -t726 - t791;
t690 = (t549 + t726) * t641;
t550 = -t725 - t789;
t689 = (t550 + t725) * t643;
t685 = t548 * t633 - t712;
t684 = t549 * t635 - t710;
t683 = t550 * t637 - t708;
t682 = -0.2e1 * qJ(2,1) * t800 - 0.2e1 * t589 * t737;
t681 = 0.2e1 * qJ(2,2) * t801 + 0.2e1 * t586 * t738;
t680 = 0.2e1 * qJ(2,3) * t802 + 0.2e1 * t583 * t739;
t676 = t524 * t775 + t525 * t772 + t526 * t769;
t656 = pkin(1) ^ 2;
t574 = t644 * t601 + t638 * t604;
t573 = t642 * t600 + t636 * t603;
t572 = t599 * t640 + t602 * t634;
t571 = -t638 * t601 + t644 * t604;
t570 = -t636 * t600 + t642 * t603;
t569 = -t599 * t634 + t602 * t640;
t568 = t594 * t644 + t613 * t638;
t567 = -t593 * t642 + t613 * t636;
t566 = -t592 * t640 + t613 * t634;
t565 = -t594 * t638 + t613 * t644;
t564 = t593 * t636 + t613 * t642;
t563 = t592 * t634 + t613 * t640;
t556 = t565 * t601 + t568 * t604;
t555 = t564 * t600 + t567 * t603;
t554 = t563 * t599 + t566 * t602;
t553 = t565 * t604 - t568 * t601;
t552 = t564 * t603 - t567 * t600;
t551 = t563 * t602 - t566 * t599;
t517 = t677 - 0.2e1 * t806;
t516 = t677 - t806;
t515 = t678 - 0.2e1 * t807;
t514 = t678 - t807;
t513 = t679 - 0.2e1 * t808;
t512 = t679 - t808;
t508 = (t651 + t656) * t526 - pkin(1) * t677;
t507 = (t650 + t656) * t525 - pkin(1) * t678;
t506 = (t649 + t656) * t524 - pkin(1) * t679;
t1 = [t569 * t802 + t570 * t801 - t571 * t800, 0, 0, -(t517 * t571 + t526 * t553) * t587 + (t515 * t570 + t525 * t552) * t584 + (t513 * t569 + t524 * t551) * t581, -t551 * t792 - t552 * t790 + t553 * t788 + t569 * t680 + t570 * t681 + t571 * t682, -(t571 * t508 + t553 * t516) * t587 + (t570 * t507 + t552 * t514) * t584 + (t569 * t506 + t551 * t512) * t581 + (-t551 * t557 + t569 * t722) * t803 + (-t552 * t558 + t570 * t721) * t804 - (-t553 * t559 + t571 * t720) * t805, t569 * t742 + t570 * t741 - t571 * t740 + (t569 * t733 + t570 * t732 + t571 * t731) * t723, -0.2e1 * t569 * t719 - 0.2e1 * t570 * t718 + 0.2e1 * t571 * t717 + 0.2e1 * (t569 * t716 + t570 * t715 + t571 * t714) * t759, (t569 * t701 + t570 * t700 - t571 * t699) * t766, 0, 0, -(t683 * t553 + t753 * t571) * t587 + (t684 * t552 + t755 * t570) * t584 + (t685 * t551 + t757 * t569) * t581, -(t553 * t689 + t754 * t571) * t587 + (t552 * t690 + t756 * t570) * t584 + (t551 * t691 + t758 * t569) * t581, 0; t572 * t802 + t573 * t801 - t574 * t800, 0, 0, -(t517 * t574 + t526 * t556) * t587 + (t515 * t573 + t525 * t555) * t584 + (t513 * t572 + t524 * t554) * t581, -t554 * t792 - t555 * t790 + t556 * t788 + t572 * t680 + t573 * t681 + t574 * t682, -(t574 * t508 + t556 * t516) * t587 + (t573 * t507 + t555 * t514) * t584 + (t572 * t506 + t554 * t512) * t581 + (-t554 * t557 + t572 * t722) * t803 + (-t555 * t558 + t573 * t721) * t804 - (-t556 * t559 + t574 * t720) * t805, t572 * t742 + t573 * t741 - t574 * t740 + (t572 * t733 + t573 * t732 + t574 * t731) * t723, -0.2e1 * t572 * t719 - 0.2e1 * t573 * t718 + 0.2e1 * t574 * t717 + 0.2e1 * (t572 * t716 + t573 * t715 + t574 * t714) * t759, (t572 * t701 + t573 * t700 - t574 * t699) * t766, 0, 0, -(t683 * t556 + t753 * t574) * t587 + (t684 * t555 + t755 * t573) * t584 + (t685 * t554 + t757 * t572) * t581, -(t556 * t689 + t754 * t574) * t587 + (t555 * t690 + t756 * t573) * t584 + (t554 * t691 + t758 * t572) * t581, 0; 0, 0, 0, t676, -t614 * t736 - t618 * t735 - t622 * t734, (t516 - t748) * t769 + (t514 - t746) * t772 + (t512 - t744) * t775, (-t734 - t735 - t736) * t653, (-t557 * t730 - t558 * t729 - t559 * t728) * t653, -t676 * t653, (t524 + t525 + t526) * t653, (t767 + t770 + t773) * t653 / t652 * t629, t548 * t639 + t549 * t641 + t550 * t643 + (-t626 * t773 - t627 * t770 - t628 * t767) * t766 + (t705 * t775 + t706 * t772 + t707 * t769) * t653, t548 * t814 + t549 * t813 + t550 * t812 + (t768 + t771 + t774) * t766 + (-t705 - t706 - t707) * t653, 0;];
tau_reg  = t1;
