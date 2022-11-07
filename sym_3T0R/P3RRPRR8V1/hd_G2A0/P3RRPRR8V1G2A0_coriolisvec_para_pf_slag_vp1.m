% Calculate vector of centrifugal and Coriolis load for parallel robot
% P3RRPRR8V1G2A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4,theta3]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2022-11-04 17:05
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(5,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-04 17:04:15
% EndTime: 2022-11-04 17:04:17
% DurationCPUTime: 2.74s
% Computational Cost: add. (13494->288), mult. (22866->473), div. (2937->12), fcn. (18921->44), ass. (0->241)
t721 = pkin(4) + qJ(3,3);
t728 = sin(qJ(1,3));
t734 = cos(qJ(1,3));
t717 = cos(pkin(5));
t825 = t717 * pkin(2);
t667 = pkin(1) + t825;
t733 = cos(qJ(2,3));
t716 = sin(pkin(5));
t727 = sin(qJ(2,3));
t809 = t716 * t727;
t759 = pkin(2) * t809 - t667 * t733;
t624 = -t734 * t721 - t759 * t728;
t829 = pkin(2) * t716;
t636 = t727 * t667 + t733 * t829;
t724 = legFrame(3,2);
t686 = sin(t724);
t689 = cos(t724);
t602 = t624 * t689 + t686 * t636;
t603 = -t624 * t686 + t689 * t636;
t703 = qJ(2,3) + pkin(5);
t682 = cos(t703);
t696 = t733 * pkin(1);
t850 = pkin(2) * t682 + t696;
t627 = t728 * t721 + t734 * t850;
t707 = 0.1e1 / t721;
t744 = xDP(3);
t745 = xDP(2);
t746 = xDP(1);
t590 = (t602 * t746 + t603 * t745 + t627 * t744) * t707;
t863 = 0.2e1 * t590;
t722 = pkin(4) + qJ(3,2);
t730 = sin(qJ(1,2));
t736 = cos(qJ(1,2));
t735 = cos(qJ(2,2));
t729 = sin(qJ(2,2));
t808 = t716 * t729;
t758 = pkin(2) * t808 - t667 * t735;
t625 = -t736 * t722 - t758 * t730;
t637 = t729 * t667 + t735 * t829;
t725 = legFrame(2,2);
t687 = sin(t725);
t690 = cos(t725);
t604 = t625 * t690 + t687 * t637;
t605 = -t625 * t687 + t690 * t637;
t705 = qJ(2,2) + pkin(5);
t684 = cos(t705);
t697 = t735 * pkin(1);
t849 = pkin(2) * t684 + t697;
t628 = t730 * t722 + t736 * t849;
t708 = 0.1e1 / t722;
t591 = (t604 * t746 + t605 * t745 + t628 * t744) * t708;
t862 = 0.2e1 * t591;
t723 = pkin(4) + qJ(3,1);
t732 = sin(qJ(1,1));
t738 = cos(qJ(1,1));
t737 = cos(qJ(2,1));
t731 = sin(qJ(2,1));
t807 = t716 * t731;
t757 = pkin(2) * t807 - t667 * t737;
t626 = -t738 * t723 - t757 * t732;
t638 = t731 * t667 + t737 * t829;
t726 = legFrame(1,2);
t688 = sin(t726);
t691 = cos(t726);
t606 = t626 * t691 + t688 * t638;
t607 = -t626 * t688 + t691 * t638;
t706 = qJ(2,1) + pkin(5);
t685 = cos(t706);
t698 = t737 * pkin(1);
t848 = pkin(2) * t685 + t698;
t629 = t732 * t723 + t738 * t848;
t709 = 0.1e1 / t723;
t592 = (t606 * t746 + t607 * t745 + t629 * t744) * t709;
t861 = 0.2e1 * t592;
t797 = 0.2e1 * pkin(1);
t639 = 0.1e1 / t850;
t617 = (t686 * t746 + t689 * t745) * t639;
t614 = t617 ^ 2;
t754 = pkin(2) ^ 2;
t755 = pkin(1) ^ 2;
t798 = -t754 - t755;
t649 = t825 * t797 - t798;
t859 = t614 * t649;
t640 = 0.1e1 / t849;
t618 = (t687 * t746 + t690 * t745) * t640;
t615 = t618 ^ 2;
t858 = t615 * t649;
t641 = 0.1e1 / t848;
t619 = (t688 * t746 + t691 * t745) * t641;
t616 = t619 ^ 2;
t857 = t616 * t649;
t676 = sin(t703);
t693 = t727 * pkin(1);
t856 = t676 * rSges(3,1) + t682 * rSges(3,2) + t693;
t855 = pkin(2) * t676 + t693;
t678 = sin(t705);
t694 = t729 * pkin(1);
t854 = t678 * rSges(3,1) + t684 * rSges(3,2) + t694;
t853 = pkin(2) * t678 + t694;
t679 = sin(t706);
t695 = t731 * pkin(1);
t852 = t679 * rSges(3,1) + t685 * rSges(3,2) + t695;
t851 = pkin(2) * t679 + t695;
t847 = 2 * m(3);
t846 = pkin(1) * m(3);
t845 = m(2) * rSges(2,3);
t844 = m(3) * rSges(3,2);
t751 = rSges(2,2) ^ 2;
t753 = rSges(2,1) ^ 2;
t642 = t755 * m(3) + (-t751 + t753) * m(2) + Icges(2,2) - Icges(2,1);
t843 = -t642 / 0.2e1;
t750 = rSges(3,2) ^ 2;
t752 = rSges(3,1) ^ 2;
t651 = (-t750 + t752) * m(3) - Icges(3,1) + Icges(3,2);
t842 = -t651 / 0.2e1;
t841 = t651 / 0.2e1;
t671 = -rSges(3,1) * t844 + Icges(3,4);
t840 = -t671 / 0.2e1;
t673 = m(2) * rSges(2,1) * rSges(2,2) - Icges(2,4);
t839 = t673 / 0.2e1;
t718 = (rSges(3,3) + qJ(3,3));
t838 = m(3) * t718;
t719 = (rSges(3,3) + qJ(3,2));
t837 = m(3) * t719;
t720 = (rSges(3,3) + qJ(3,1));
t836 = m(3) * t720;
t828 = t590 * m(3);
t827 = t591 * m(3);
t826 = t592 * m(3);
t747 = 0.2e1 * qJ(2,3);
t710 = sin(t747);
t824 = t642 * t710;
t748 = 0.2e1 * qJ(2,2);
t711 = sin(t748);
t823 = t642 * t711;
t749 = 0.2e1 * qJ(2,1);
t712 = sin(t749);
t822 = t642 * t712;
t821 = t667 * t689;
t820 = t667 * t690;
t819 = t667 * t691;
t668 = 0.2e1 * t703;
t657 = cos(t668);
t818 = t671 * t657;
t669 = 0.2e1 * t705;
t658 = cos(t669);
t817 = t671 * t658;
t670 = 0.2e1 * t706;
t659 = cos(t670);
t816 = t671 * t659;
t713 = cos(t747);
t815 = t673 * t713;
t714 = cos(t748);
t814 = t673 * t714;
t715 = cos(t749);
t813 = t673 * t715;
t812 = t686 * t667;
t811 = t687 * t667;
t810 = t688 * t667;
t806 = -rSges(2,1) * t845 + Icges(2,5);
t672 = rSges(2,2) * t845 - Icges(2,6);
t792 = t686 * t829;
t795 = t689 * t829;
t608 = (-t728 * t812 + t795) * t733 + (t728 * t792 + t821) * t727;
t611 = (t728 * t821 + t792) * t733 + (-t728 * t795 + t812) * t727;
t633 = 0.1e1 / t759;
t587 = (t734 * t744 - (t608 * t745 + t611 * t746) * t633) * t707;
t788 = pkin(2) * t797;
t701 = pkin(5) + t747;
t680 = cos(t701);
t802 = -t680 - t717;
t569 = (-t859 + ((-t754 * t657 - t755 * t713 + t802 * t788 + t798) * t587 / 0.2e1 + t850 * t863 + (-t587 * t721 + 0.2e1 * t617 * t855) * t721) * t587) * t707;
t581 = (0.1e1 / (-t696 + (-t717 * t733 + t809) * pkin(2)) * t859 + (t759 * t587 + t863) * t587) * t707;
t660 = rSges(3,2) * t838 - Icges(3,6);
t663 = rSges(3,1) * t838 - Icges(3,5);
t765 = pkin(1) * t838 - t806;
t593 = -(t660 * t716 - t663 * t717 - t765) * t727 + (t660 * t717 + t663 * t716 + t672) * t733;
t630 = t682 * rSges(3,1) - t676 * rSges(3,2) + t696;
t654 = sin(t668);
t674 = sin(t701);
t666 = pkin(1) * t716 * t844;
t799 = t751 + t753;
t756 = -(rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - ((2 * rSges(2,3) ^ 2) + t799) * m(2) / 0.2e1 + t666 - Icges(3,2) / 0.2e1 - Icges(2,2) / 0.2e1 - Icges(3,1) / 0.2e1 - Icges(2,1) / 0.2e1 - Icges(1,3);
t771 = rSges(3,1) * t674 + rSges(3,2) * t680;
t772 = -t750 / 0.2e1 - t752 / 0.2e1 - t755 / 0.2e1;
t787 = t614 * t639 * t855;
t789 = m(3) * t797;
t805 = (t657 * t842 + t713 * t843 + t756) * t581 - t593 * t787 + 0.2e1 * (t654 * t840 + t710 * t839) * t581 + (t630 * t569 + (-(t718 ^ 2) + t772 + (t802 * rSges(3,1) + rSges(3,2) * t674) * pkin(1)) * t581) * m(3) + (t660 * t676 - t663 * t682 + t672 * t727 - t765 * t733) * t614 + (0.2e1 * t718 * t828 + (-t651 * t654 - t771 * t789 - 0.2e1 * t815 + 0.2e1 * t818 - t824) * t617) * t587;
t791 = t687 * t829;
t794 = t690 * t829;
t609 = (-t730 * t811 + t794) * t735 + (t730 * t791 + t820) * t729;
t612 = (t730 * t820 + t791) * t735 + (-t730 * t794 + t811) * t729;
t634 = 0.1e1 / t758;
t588 = (t736 * t744 - (t609 * t745 + t612 * t746) * t634) * t708;
t704 = t748 + pkin(5);
t683 = cos(t704);
t800 = -t683 - t717;
t570 = (-t858 + ((-t754 * t658 - t755 * t714 + t800 * t788 + t798) * t588 / 0.2e1 + t849 * t862 + (-t588 * t722 + 0.2e1 * t618 * t853) * t722) * t588) * t708;
t582 = (-0.1e1 / (t697 + (t717 * t735 - t808) * pkin(2)) * t858 + (t758 * t588 + t862) * t588) * t708;
t661 = rSges(3,2) * t837 - Icges(3,6);
t664 = rSges(3,1) * t837 - Icges(3,5);
t764 = pkin(1) * t837 - t806;
t594 = -(t661 * t716 - t664 * t717 - t764) * t729 + (t661 * t717 + t664 * t716 + t672) * t735;
t631 = t684 * rSges(3,1) - t678 * rSges(3,2) + t697;
t655 = sin(t669);
t677 = sin(t704);
t768 = rSges(3,1) * t677 + rSges(3,2) * t683;
t786 = t615 * t640 * t853;
t804 = (t658 * t842 + t714 * t843 + t756) * t582 - t594 * t786 + 0.2e1 * (t655 * t840 + t711 * t839) * t582 + (t631 * t570 + (-(t719 ^ 2) + t772 + (t800 * rSges(3,1) + rSges(3,2) * t677) * pkin(1)) * t582) * m(3) + (t661 * t678 - t664 * t684 + t672 * t729 - t764 * t735) * t615 + (0.2e1 * t719 * t827 + (-t651 * t655 - t768 * t789 - 0.2e1 * t814 + 0.2e1 * t817 - t823) * t618) * t588;
t790 = t688 * t829;
t793 = t691 * t829;
t610 = (-t732 * t810 + t793) * t737 + (t732 * t790 + t819) * t731;
t613 = (t732 * t819 + t790) * t737 + (-t732 * t793 + t810) * t731;
t635 = 0.1e1 / t757;
t589 = (t738 * t744 - (t610 * t745 + t613 * t746) * t635) * t709;
t702 = pkin(5) + t749;
t681 = cos(t702);
t801 = -t681 - t717;
t571 = (-t857 + ((-t754 * t659 - t755 * t715 + t801 * t788 + t798) * t589 / 0.2e1 + t848 * t861 + (-t589 * t723 + 0.2e1 * t619 * t851) * t723) * t589) * t709;
t583 = (-0.1e1 / (t698 + (t717 * t737 - t807) * pkin(2)) * t857 + (t757 * t589 + t861) * t589) * t709;
t662 = rSges(3,2) * t836 - Icges(3,6);
t665 = rSges(3,1) * t836 - Icges(3,5);
t763 = pkin(1) * t836 - t806;
t595 = -(t662 * t716 - t665 * t717 - t763) * t731 + (t662 * t717 + t665 * t716 + t672) * t737;
t632 = t685 * rSges(3,1) - t679 * rSges(3,2) + t698;
t656 = sin(t670);
t675 = sin(t702);
t770 = rSges(3,1) * t675 + rSges(3,2) * t681;
t785 = t616 * t641 * t851;
t803 = (t659 * t842 + t715 * t843 + t756) * t583 - t595 * t785 + 0.2e1 * (t656 * t840 + t712 * t839) * t583 + (t632 * t571 + (-(t720 ^ 2) + t772 + (t801 * rSges(3,1) + rSges(3,2) * t675) * pkin(1)) * t583) * m(3) + (t662 * t679 - t665 * t685 + t672 * t731 - t763 * t737) * t616 + (0.2e1 * t720 * t826 + (-t651 * t656 - t770 * t789 - 0.2e1 * t813 + 0.2e1 * t816 - t822) * t619) * t589;
t784 = t805 * t633;
t783 = t804 * t634;
t782 = t803 * t635;
t781 = (-t587 * t718 / 0.2e1 + t856 * t617) * t587 * t847 + (t581 * t630 - t569) * m(3);
t780 = (-t588 * t719 / 0.2e1 + t854 * t618) * t588 * t847 + (t582 * t631 - t570) * m(3);
t779 = (-t589 * t720 / 0.2e1 + t852 * t619) * t589 * t847 + (t583 * t632 - t571) * m(3);
t620 = 0.2e1 * t666 - t799 * m(2) - Icges(2,3) - Icges(3,3) + (-0.2e1 * pkin(1) * rSges(3,1) * t717 - t750 - t752 - t755) * m(3);
t762 = t639 * ((-0.2e1 * t856 * t828 + (t654 * t841 - t818 + t824 / 0.2e1 + t815 + t771 * t846) * t587) * t587 + t593 * t581 - t620 * t787);
t761 = t640 * ((-0.2e1 * t854 * t827 + (t655 * t841 - t817 + t823 / 0.2e1 + t814 + t768 * t846) * t588) * t588 + t594 * t582 - t620 * t786);
t760 = t641 * ((-0.2e1 * t852 * t826 + (t656 * t841 - t816 + t822 / 0.2e1 + t813 + t770 * t846) * t589) * t589 + t595 * t583 - t620 * t785);
t1 = [t688 * t760 + t687 * t761 + t686 * t762 + (t779 * t606 - t613 * t782) * t709 + (t780 * t604 - t612 * t783) * t708 + (t781 * t602 - t611 * t784) * t707; t691 * t760 + t690 * t761 + t689 * t762 + (t779 * t607 - t610 * t782) * t709 + (t780 * t605 - t609 * t783) * t708 + (t781 * t603 - t608 * t784) * t707; (t779 * t629 + t803 * t738) * t709 + (t780 * t628 + t804 * t736) * t708 + (t781 * t627 + t805 * t734) * t707;];
taucX  = t1;
