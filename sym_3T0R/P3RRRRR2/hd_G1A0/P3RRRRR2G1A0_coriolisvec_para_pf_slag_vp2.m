% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRRRR2G1A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:05
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRRRR2G1A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G1A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRRRR2G1A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G1A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G1A0_coriolisvec_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR2G1A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR2G1A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRRRR2G1A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G1A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G1A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:04:18
% EndTime: 2020-03-09 21:04:20
% DurationCPUTime: 2.39s
% Computational Cost: add. (12930->246), mult. (9339->445), div. (7065->18), fcn. (9309->42), ass. (0->230)
t849 = 2 * pkin(1);
t685 = qJ(1,3) + legFrame(3,3);
t676 = qJ(2,3) + t685;
t670 = qJ(3,3) + t676;
t671 = -qJ(3,3) + t676;
t842 = -2 * pkin(1);
t634 = sin(t685) * t842 + (-sin(t671) - sin(t670)) * pkin(2);
t652 = 0.1e1 / (sin(qJ(2,3) + qJ(3,3)) + sin(qJ(2,3) - qJ(3,3)));
t722 = xDP(2);
t726 = 0.1e1 / pkin(2);
t729 = 1 / pkin(1);
t793 = t726 * t729;
t773 = t722 * t793;
t619 = t634 * t652 * t773;
t637 = cos(t685) * t842 + (-cos(t671) - cos(t670)) * pkin(2);
t723 = xDP(1);
t772 = t723 * t793;
t622 = t637 * t652 * t772;
t721 = xDP(3);
t715 = cos(qJ(3,3));
t709 = sin(qJ(3,3));
t710 = sin(qJ(2,3));
t688 = 0.1e1 / t710;
t819 = t688 * t729;
t785 = t709 * t819;
t693 = 0.1e1 / t715 ^ 2;
t815 = t693 * t726;
t716 = cos(qJ(2,3));
t835 = pkin(1) * t716;
t744 = (pkin(2) * t715 + t835) * t785 * t815;
t738 = t721 * t744;
t616 = t619 + t622 - t738;
t664 = sin(t676);
t667 = cos(t676);
t692 = 0.1e1 / t715;
t763 = t692 * t785;
t794 = t723 * t729;
t795 = t722 * t729;
t628 = t721 * t763 + (t664 * t795 + t667 * t794) * t688;
t613 = t628 + t616;
t725 = pkin(2) ^ 2;
t691 = t715 ^ 2;
t730 = t715 * t691;
t604 = t725 * t613 * t730;
t607 = t622 / 0.2e1 + t619 / 0.2e1 - t738 / 0.2e1 + t628;
t694 = 0.1e1 / t730;
t728 = pkin(1) ^ 2;
t796 = t721 * t726;
t782 = t692 * t796;
t760 = t716 * t782;
t799 = t709 * t710;
t776 = t721 * t799;
t832 = pkin(2) * t691;
t789 = t716 * t832;
t583 = (((-pkin(1) * t613 * t799 + t692 * t721) * t715 + pkin(1) * t760) * t694 * t796 + ((t604 + t607 * t789 * t849 + (-pkin(1) * t692 * t776 + t628 * t728) * t715) * t628 + (t604 + (t613 * t789 - t776) * pkin(1)) * t616) * t815) * t819;
t703 = t721 ^ 2;
t810 = t703 * t726;
t822 = t613 * t715;
t589 = ((-t628 * t715 * t835 - t613 * t832) * t692 * t628 - pkin(2) * t616 * t822 - t694 * t810) * t819;
t682 = mrSges(3,1) * t835;
t679 = mrSges(3,2) * t709 - mrSges(2,1);
t708 = mrSges(2,2) - mrSges(3,3);
t741 = (t679 * t716 + t708 * t710) * pkin(1);
t707 = Ifges(3,2) - Ifges(3,1);
t808 = t707 * t691;
t829 = -Ifges(3,1) - Ifges(2,3);
t754 = -t808 + t829;
t828 = Ifges(3,4) * t709;
t631 = -(t682 + 0.2e1 * t828) * t715 + t741 + t754;
t838 = pkin(1) * t710;
t640 = -t715 * (-mrSges(3,2) * t838 + Ifges(3,6)) + t709 * (mrSges(3,1) * t838 - Ifges(3,5));
t751 = -0.2e1 * t613 * t782;
t770 = 0.4e1 * Ifges(3,4) * t796;
t809 = t703 / pkin(2) ^ 2;
t779 = t709 * t809;
t786 = Ifges(3,5) * t809;
t805 = t707 * t709;
t735 = t770 * t822 + Ifges(3,4) * t751 + (t693 * t786 + t751 * t805) * t715 - Ifges(3,6) * t693 * t779;
t745 = -((m(2) + m(3)) * t728) - Ifges(1,3) + t829;
t748 = t613 * t760;
t759 = t694 * t779;
t802 = t708 * t716;
t771 = t809 / 0.2e1;
t825 = (t693 * t771 + (t628 + t616 / 0.2e1) * t616) * t710;
t841 = -0.2e1 * t715;
t848 = t735 + ((-mrSges(3,1) * t825 - mrSges(3,2) * t748) * t715 + (-mrSges(3,1) * t748 + mrSges(3,2) * t825) * t709 + (-mrSges(2,1) * t710 - t802) * t616 * t607) * t849 + (-t808 + (t682 + t828) * t841 + 0.2e1 * t741 + t745) * t589 + t631 * t583 - t640 * t759;
t686 = qJ(1,2) + legFrame(2,3);
t677 = qJ(2,2) + t686;
t672 = qJ(3,2) + t677;
t673 = -qJ(3,2) + t677;
t635 = sin(t686) * t842 + (-sin(t673) - sin(t672)) * pkin(2);
t653 = 0.1e1 / (sin(qJ(2,2) + qJ(3,2)) + sin(qJ(2,2) - qJ(3,2)));
t620 = t635 * t653 * t773;
t638 = cos(t686) * t842 + (-cos(t673) - cos(t672)) * pkin(2);
t623 = t638 * t653 * t772;
t717 = cos(qJ(3,2));
t711 = sin(qJ(3,2));
t712 = sin(qJ(2,2));
t689 = 0.1e1 / t712;
t818 = t689 * t729;
t784 = t711 * t818;
t697 = 0.1e1 / t717 ^ 2;
t813 = t697 * t726;
t718 = cos(qJ(2,2));
t834 = pkin(1) * t718;
t743 = (pkin(2) * t717 + t834) * t784 * t813;
t737 = t721 * t743;
t617 = t620 + t623 - t737;
t665 = sin(t677);
t668 = cos(t677);
t696 = 0.1e1 / t717;
t762 = t696 * t784;
t629 = t721 * t762 + (t665 * t795 + t668 * t794) * t689;
t614 = t629 + t617;
t695 = t717 ^ 2;
t731 = t717 * t695;
t605 = t725 * t614 * t731;
t608 = t623 / 0.2e1 + t620 / 0.2e1 - t737 / 0.2e1 + t629;
t698 = 0.1e1 / t731;
t781 = t696 * t796;
t758 = t718 * t781;
t798 = t711 * t712;
t775 = t721 * t798;
t831 = pkin(2) * t695;
t788 = t718 * t831;
t584 = (((-pkin(1) * t614 * t798 + t696 * t721) * t717 + pkin(1) * t758) * t698 * t796 + ((t605 + t608 * t788 * t849 + (-pkin(1) * t696 * t775 + t629 * t728) * t717) * t629 + (t605 + (t614 * t788 - t775) * pkin(1)) * t617) * t813) * t818;
t821 = t614 * t717;
t590 = ((-t629 * t717 * t834 - t614 * t831) * t696 * t629 - pkin(2) * t617 * t821 - t698 * t810) * t818;
t683 = mrSges(3,1) * t834;
t680 = mrSges(3,2) * t711 - mrSges(2,1);
t740 = (t680 * t718 + t708 * t712) * pkin(1);
t807 = t707 * t695;
t753 = -t807 + t829;
t827 = Ifges(3,4) * t711;
t632 = -(t683 + 0.2e1 * t827) * t717 + t740 + t753;
t837 = pkin(1) * t712;
t641 = -t717 * (-mrSges(3,2) * t837 + Ifges(3,6)) + t711 * (mrSges(3,1) * t837 - Ifges(3,5));
t750 = -0.2e1 * t614 * t781;
t778 = t711 * t809;
t804 = t707 * t711;
t734 = t770 * t821 + Ifges(3,4) * t750 + (t697 * t786 + t750 * t804) * t717 - Ifges(3,6) * t697 * t778;
t747 = t614 * t758;
t757 = t698 * t778;
t801 = t708 * t718;
t824 = (t697 * t771 + (t629 + t617 / 0.2e1) * t617) * t712;
t840 = -0.2e1 * t717;
t847 = t734 + ((-mrSges(3,1) * t824 - mrSges(3,2) * t747) * t717 + (-mrSges(3,1) * t747 + mrSges(3,2) * t824) * t711 + (-mrSges(2,1) * t712 - t801) * t617 * t608) * t849 + (-t807 + (t683 + t827) * t840 + 0.2e1 * t740 + t745) * t590 + t632 * t584 - t641 * t757;
t687 = qJ(1,1) + legFrame(1,3);
t678 = qJ(2,1) + t687;
t674 = qJ(3,1) + t678;
t675 = -qJ(3,1) + t678;
t636 = sin(t687) * t842 + (-sin(t675) - sin(t674)) * pkin(2);
t654 = 0.1e1 / (sin(qJ(2,1) + qJ(3,1)) + sin(qJ(2,1) - qJ(3,1)));
t621 = t636 * t654 * t773;
t639 = cos(t687) * t842 + (-cos(t675) - cos(t674)) * pkin(2);
t624 = t639 * t654 * t772;
t719 = cos(qJ(3,1));
t713 = sin(qJ(3,1));
t714 = sin(qJ(2,1));
t690 = 0.1e1 / t714;
t817 = t690 * t729;
t783 = t713 * t817;
t701 = 0.1e1 / t719 ^ 2;
t811 = t701 * t726;
t720 = cos(qJ(2,1));
t833 = pkin(1) * t720;
t742 = (pkin(2) * t719 + t833) * t783 * t811;
t736 = t721 * t742;
t618 = t621 + t624 - t736;
t666 = sin(t678);
t669 = cos(t678);
t700 = 0.1e1 / t719;
t761 = t700 * t783;
t630 = t721 * t761 + (t666 * t795 + t669 * t794) * t690;
t615 = t618 + t630;
t699 = t719 ^ 2;
t732 = t719 * t699;
t606 = t725 * t615 * t732;
t609 = t624 / 0.2e1 + t621 / 0.2e1 - t736 / 0.2e1 + t630;
t702 = 0.1e1 / t732;
t780 = t700 * t796;
t756 = t720 * t780;
t797 = t713 * t714;
t774 = t721 * t797;
t830 = pkin(2) * t699;
t787 = t720 * t830;
t585 = (((-pkin(1) * t615 * t797 + t700 * t721) * t719 + pkin(1) * t756) * t702 * t796 + ((t606 + t609 * t787 * t849 + (-pkin(1) * t700 * t774 + t630 * t728) * t719) * t630 + (t606 + (t615 * t787 - t774) * pkin(1)) * t618) * t811) * t817;
t820 = t615 * t719;
t591 = ((-t630 * t719 * t833 - t615 * t830) * t700 * t630 - pkin(2) * t618 * t820 - t702 * t810) * t817;
t684 = mrSges(3,1) * t833;
t681 = mrSges(3,2) * t713 - mrSges(2,1);
t739 = (t681 * t720 + t708 * t714) * pkin(1);
t806 = t707 * t699;
t752 = -t806 + t829;
t826 = Ifges(3,4) * t713;
t633 = -(t684 + 0.2e1 * t826) * t719 + t739 + t752;
t836 = pkin(1) * t714;
t642 = -t719 * (-mrSges(3,2) * t836 + Ifges(3,6)) + t713 * (mrSges(3,1) * t836 - Ifges(3,5));
t749 = -0.2e1 * t615 * t780;
t777 = t713 * t809;
t803 = t707 * t713;
t733 = t770 * t820 + Ifges(3,4) * t749 + (t701 * t786 + t749 * t803) * t719 - Ifges(3,6) * t701 * t777;
t746 = t615 * t756;
t755 = t702 * t777;
t800 = t708 * t720;
t823 = (t701 * t771 + (t630 + t618 / 0.2e1) * t618) * t714;
t839 = -0.2e1 * t719;
t846 = t733 + ((-mrSges(3,1) * t823 - mrSges(3,2) * t746) * t719 + (-mrSges(3,1) * t746 + mrSges(3,2) * t823) * t713 + (-mrSges(2,1) * t714 - t800) * t618 * t609) * t849 + (-t806 + (t684 + t826) * t839 + 0.2e1 * t739 + t745) * t591 + t633 * t585 - t642 * t755;
t625 = t628 ^ 2;
t655 = -Ifges(3,5) * t709 - Ifges(3,6) * t715;
t845 = t631 * t589 + (t828 * t841 + t754) * t583 - t655 * t759 + (t802 + (mrSges(3,1) * t715 - t679) * t710) * t625 * pkin(1) + t735;
t626 = t629 ^ 2;
t656 = -Ifges(3,5) * t711 - Ifges(3,6) * t717;
t844 = t632 * t590 + (t827 * t840 + t753) * t584 - t656 * t757 + (t801 + (mrSges(3,1) * t717 - t680) * t712) * t626 * pkin(1) + t734;
t627 = t630 ^ 2;
t657 = -Ifges(3,5) * t713 - Ifges(3,6) * t719;
t843 = t633 * t591 + (t826 * t839 + t752) * t585 - t657 * t755 + (t800 + (mrSges(3,1) * t719 - t681) * t714) * t627 * pkin(1) + t733;
t792 = t625 * t835;
t791 = t626 * t834;
t790 = t627 * t833;
t769 = t848 * t688;
t768 = t847 * t689;
t767 = t846 * t690;
t766 = t845 * t652;
t765 = t844 * t653;
t764 = t843 * t654;
t612 = t615 ^ 2;
t611 = t614 ^ 2;
t610 = t613 ^ 2;
t1 = [(t669 * t767 + t668 * t768 + t667 * t769 + (t637 * t766 + t638 * t765 + t639 * t764) * t726) * t729; (t666 * t767 + t665 * t768 + t664 * t769 + (t634 * t766 + t635 * t765 + t636 * t764) * t726) * t729; -t843 * t742 - t844 * t743 - t845 * t744 + t846 * t761 + t847 * t762 + t848 * t763 + ((Ifges(3,3) * t759 + t655 * t583 + t640 * t589 + (mrSges(3,2) * t792 + t610 * t805) * t715 + mrSges(3,1) * t709 * t792 + (-0.2e1 * t691 + 0.1e1) * t610 * Ifges(3,4)) * t692 + (Ifges(3,3) * t757 + t656 * t584 + t641 * t590 + (mrSges(3,2) * t791 + t611 * t804) * t717 + mrSges(3,1) * t711 * t791 + (-0.2e1 * t695 + 0.1e1) * t611 * Ifges(3,4)) * t696 + (Ifges(3,3) * t755 + t657 * t585 + t642 * t591 + (mrSges(3,2) * t790 + t612 * t803) * t719 + mrSges(3,1) * t713 * t790 + (-0.2e1 * t699 + 0.1e1) * t612 * Ifges(3,4)) * t700) * t726;];
taucX  = t1;
