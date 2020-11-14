% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRRR8V2G3A0
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2020-08-06 18:05
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G3A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:04:36
% EndTime: 2020-08-06 18:04:41
% DurationCPUTime: 5.44s
% Computational Cost: add. (33657->292), mult. (75984->556), div. (4716->7), fcn. (74370->22), ass. (0->202)
t733 = sin(qJ(2,1));
t739 = cos(qJ(2,1));
t745 = pkin(7) + pkin(6);
t679 = pkin(2) * t733 - t739 * t745;
t721 = sin(pkin(4));
t723 = cos(pkin(4));
t732 = sin(qJ(3,1));
t791 = t732 * t723;
t664 = pkin(3) * t791 + t721 * t679;
t738 = cos(qJ(3,1));
t810 = t721 * t733;
t719 = t738 ^ 2;
t836 = pkin(3) * t719;
t643 = 0.1e1 / (pkin(2) * t791 + t664 * t738 + t810 * t836);
t731 = sin(qJ(2,2));
t737 = cos(qJ(2,2));
t678 = pkin(2) * t731 - t737 * t745;
t730 = sin(qJ(3,2));
t793 = t730 * t723;
t663 = pkin(3) * t793 + t721 * t678;
t736 = cos(qJ(3,2));
t812 = t721 * t731;
t718 = t736 ^ 2;
t837 = pkin(3) * t718;
t642 = 0.1e1 / (pkin(2) * t793 + t663 * t736 + t812 * t837);
t729 = sin(qJ(2,3));
t735 = cos(qJ(2,3));
t677 = pkin(2) * t729 - t735 * t745;
t728 = sin(qJ(3,3));
t795 = t728 * t723;
t662 = pkin(3) * t795 + t721 * t677;
t734 = cos(qJ(3,3));
t814 = t721 * t729;
t717 = t734 ^ 2;
t838 = pkin(3) * t717;
t641 = 0.1e1 / (pkin(2) * t795 + t662 * t734 + t814 * t838);
t680 = pkin(2) * t735 + t729 * t745;
t720 = sin(pkin(8));
t722 = cos(pkin(8));
t802 = t723 * t735;
t806 = t722 * t723;
t835 = t734 * pkin(3);
t635 = (t720 * t729 - t722 * t802) * t835 - t680 * t806 + t677 * t720;
t805 = t723 * t729;
t671 = t720 * t735 + t722 * t805;
t809 = t721 * t734;
t650 = t728 * t671 + t722 * t809;
t668 = t720 * t805 - t722 * t735;
t647 = t728 * t668 + t720 * t809;
t742 = xDP(3);
t725 = legFrame(3,2);
t705 = sin(t725);
t708 = cos(t725);
t743 = xDP(2);
t744 = xDP(1);
t757 = t705 * t743 - t708 * t744;
t623 = (t647 * t742 + t650 * t757) * t641;
t823 = t623 * t745;
t772 = t728 * t823;
t816 = t720 * t723;
t638 = (t720 * t802 + t722 * t729) * t835 + t680 * t816 + t677 * t722;
t749 = 0.1e1 / pkin(3);
t617 = (-t635 * t757 + t638 * t742) * t749 * t641;
t841 = pkin(3) * t617;
t608 = t772 - t841;
t775 = t728 * t841;
t789 = t735 * t623;
t794 = t729 * t734;
t815 = t721 * t728;
t590 = (((t723 * t617 + t721 * t789) * t838 + ((-t775 + t823) * t729 + pkin(2) * t789) * t809 + t723 * t608) * t623 + (t735 * t721 * t617 + (t717 * t723 - t794 * t815 - t723) * t623) * t841) * t641;
t665 = pkin(3) * t794 + t677;
t799 = t723 * t749;
t750 = pkin(2) ^ 2;
t698 = t745 ^ 2 + t750;
t748 = pkin(3) ^ 2;
t786 = 0.2e1 * pkin(3) * pkin(2);
t826 = (-t745 * t775 + (t717 * t748 + t734 * t786 + t698) * t623) * t623;
t593 = t641 * t799 * t826 + (-t723 * t772 + (-t665 * t815 + t723 * (pkin(2) * t734 + t838)) * t617) / (t665 * t809 + (pkin(2) + t835) * t795) * t617;
t599 = (-t734 * t826 - (pkin(2) * t617 - t608 * t734) * t841) * t641;
t620 = t623 ^ 2;
t829 = t728 * mrSges(3,1);
t686 = t734 * mrSges(3,2) + t829;
t769 = t734 * mrSges(3,1) - mrSges(3,2) * t728;
t653 = t686 * t814 - t723 * t769;
t700 = pkin(6) * mrSges(3,2) - Ifges(3,6);
t701 = mrSges(3,1) * pkin(6) - Ifges(3,5);
t674 = t700 * t734 + t728 * t701;
t724 = Ifges(3,1) - Ifges(3,2);
t740 = mrSges(3,2) * pkin(2);
t847 = -2 * Ifges(3,4);
t848 = Ifges(3,4) + t717 * t847 + (-t724 * t728 + t740) * t734;
t760 = (-Ifges(3,3) * t593 + t674 * t590 + t653 * t599 + t620 * (pkin(2) * t829 + t848)) * t749;
t741 = mrSges(3,1) * pkin(2);
t751 = -0.2e1 * pkin(6) * mrSges(3,3) + (-pkin(6) ^ 2 - t750) * m(3) - Ifges(3,1) - Ifges(2,3);
t779 = 0.2e1 * t740;
t699 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t711 = m(3) * pkin(2) + mrSges(2,1);
t820 = ((t769 + t711) * t735 + t699 * t729) * t721;
t845 = t701 / 0.2e1;
t846 = -t700 / 0.2e1;
t766 = 0.2e1 * ((t846 * t728 + t845 * t734) * t617 + (t741 * t728 + t848) * t623) * t617 + t599 * t820 - (t724 * t717 - 0.2e1 * (Ifges(3,4) * t728 + t741) * t734 + t728 * t779 + t751) * t590 - t674 * t593;
t862 = t635 * t760 + t766 * t650;
t681 = pkin(2) * t737 + t731 * t745;
t801 = t723 * t737;
t834 = t736 * pkin(3);
t636 = (t720 * t731 - t722 * t801) * t834 - t681 * t806 + t678 * t720;
t804 = t723 * t731;
t672 = t720 * t737 + t722 * t804;
t808 = t721 * t736;
t651 = t730 * t672 + t722 * t808;
t669 = t720 * t804 - t722 * t737;
t648 = t730 * t669 + t720 * t808;
t726 = legFrame(2,2);
t706 = sin(t726);
t709 = cos(t726);
t756 = t706 * t743 - t709 * t744;
t624 = (t648 * t742 + t651 * t756) * t642;
t822 = t624 * t745;
t771 = t730 * t822;
t639 = (t720 * t801 + t722 * t731) * t834 + t681 * t816 + t678 * t722;
t618 = (-t636 * t756 + t639 * t742) * t749 * t642;
t840 = pkin(3) * t618;
t609 = t771 - t840;
t774 = t730 * t840;
t788 = t737 * t624;
t792 = t731 * t736;
t813 = t721 * t730;
t591 = (((t723 * t618 + t721 * t788) * t837 + ((-t774 + t822) * t731 + pkin(2) * t788) * t808 + t723 * t609) * t624 + (t737 * t721 * t618 + (t718 * t723 - t792 * t813 - t723) * t624) * t840) * t642;
t666 = pkin(3) * t792 + t678;
t825 = (-t745 * t774 + (t718 * t748 + t736 * t786 + t698) * t624) * t624;
t594 = t642 * t799 * t825 + (-t723 * t771 + (-t666 * t813 + t723 * (pkin(2) * t736 + t837)) * t618) / (t666 * t808 + (pkin(2) + t834) * t793) * t618;
t600 = (-t736 * t825 - (pkin(2) * t618 - t609 * t736) * t840) * t642;
t621 = t624 ^ 2;
t828 = t730 * mrSges(3,1);
t687 = t736 * mrSges(3,2) + t828;
t768 = t736 * mrSges(3,1) - mrSges(3,2) * t730;
t654 = t687 * t812 - t723 * t768;
t675 = t700 * t736 + t730 * t701;
t849 = Ifges(3,4) + t718 * t847 + (-t724 * t730 + t740) * t736;
t759 = (-Ifges(3,3) * t594 + t675 * t591 + t654 * t600 + t621 * (pkin(2) * t828 + t849)) * t749;
t819 = ((t768 + t711) * t737 + t699 * t731) * t721;
t765 = 0.2e1 * ((t846 * t730 + t845 * t736) * t618 + (t741 * t730 + t849) * t624) * t618 + t600 * t819 - (t724 * t718 - 0.2e1 * (Ifges(3,4) * t730 + t741) * t736 + t730 * t779 + t751) * t591 - t675 * t594;
t861 = t636 * t759 + t765 * t651;
t682 = pkin(2) * t739 + t733 * t745;
t800 = t723 * t739;
t833 = t738 * pkin(3);
t637 = (t720 * t733 - t722 * t800) * t833 - t682 * t806 + t679 * t720;
t803 = t723 * t733;
t673 = t720 * t739 + t722 * t803;
t807 = t721 * t738;
t652 = t732 * t673 + t722 * t807;
t670 = t720 * t803 - t722 * t739;
t649 = t732 * t670 + t720 * t807;
t727 = legFrame(1,2);
t707 = sin(t727);
t710 = cos(t727);
t755 = t707 * t743 - t710 * t744;
t625 = (t649 * t742 + t652 * t755) * t643;
t821 = t625 * t745;
t770 = t732 * t821;
t640 = (t720 * t800 + t722 * t733) * t833 + t682 * t816 + t679 * t722;
t619 = (-t637 * t755 + t640 * t742) * t749 * t643;
t839 = pkin(3) * t619;
t610 = t770 - t839;
t773 = t732 * t839;
t787 = t739 * t625;
t790 = t733 * t738;
t811 = t721 * t732;
t592 = (((t723 * t619 + t721 * t787) * t836 + ((-t773 + t821) * t733 + pkin(2) * t787) * t807 + t723 * t610) * t625 + (t739 * t721 * t619 + (t719 * t723 - t790 * t811 - t723) * t625) * t839) * t643;
t667 = pkin(3) * t790 + t679;
t824 = (-t745 * t773 + (t719 * t748 + t738 * t786 + t698) * t625) * t625;
t595 = t643 * t799 * t824 + (-t723 * t770 + (-t667 * t811 + t723 * (pkin(2) * t738 + t836)) * t619) / (t667 * t807 + (pkin(2) + t833) * t791) * t619;
t601 = (-t738 * t824 - (pkin(2) * t619 - t610 * t738) * t839) * t643;
t622 = t625 ^ 2;
t827 = t732 * mrSges(3,1);
t688 = t738 * mrSges(3,2) + t827;
t767 = t738 * mrSges(3,1) - mrSges(3,2) * t732;
t655 = t688 * t810 - t723 * t767;
t676 = t700 * t738 + t732 * t701;
t850 = Ifges(3,4) + t719 * t847 + (-t724 * t732 + t740) * t738;
t758 = (-Ifges(3,3) * t595 + t676 * t592 + t655 * t601 + t622 * (pkin(2) * t827 + t850)) * t749;
t818 = ((t767 + t711) * t739 + t699 * t733) * t721;
t764 = 0.2e1 * ((t846 * t732 + t845 * t738) * t619 + (t741 * t732 + t850) * t625) * t619 + t601 * t818 - (t724 * t719 - 0.2e1 * (Ifges(3,4) * t732 + t741) * t738 + t732 * t779 + t751) * t592 - t676 * t595;
t860 = t637 * t758 + t764 * t652;
t844 = pkin(2) * t728;
t843 = pkin(2) * t730;
t842 = pkin(2) * t732;
t817 = t720 * t721;
t614 = t617 ^ 2;
t715 = -m(1) - m(2) - m(3);
t785 = -t590 * t820 + t653 * t593 + t715 * t599 + ((-t620 * t711 - t769 * (t620 + t614)) * t729 + (-0.2e1 * t617 * t686 + t623 * t699) * t789) * t721 - t614 * t723 * t686;
t615 = t618 ^ 2;
t784 = -t591 * t819 + t654 * t594 + t715 * t600 + ((-t621 * t711 - t768 * (t621 + t615)) * t731 + (-0.2e1 * t618 * t687 + t624 * t699) * t788) * t721 - t615 * t723 * t687;
t616 = t619 ^ 2;
t783 = -t592 * t818 + t655 * t595 + t715 * t601 + ((-t622 * t711 - t767 * (t622 + t616)) * t733 + (-0.2e1 * t619 * t688 + t625 * t699) * t787) * t721 - t616 * t723 * t688;
t754 = pkin(3) * t815 - t677 * t723;
t753 = pkin(3) * t813 - t678 * t723;
t752 = pkin(3) * t811 - t679 * t723;
t646 = t722 * t682 + t720 * t752;
t645 = t722 * t681 + t720 * t753;
t644 = t722 * t680 + t720 * t754;
t1 = [(t783 * (-(t670 * t710 - t707 * t810) * t836 + (t646 * t710 + t707 * t664) * t738 + (t723 * t707 + t710 * t817) * t842) + t860 * t710) * t643 + (t784 * (-(t669 * t709 - t706 * t812) * t837 + (t645 * t709 + t706 * t663) * t736 + (t723 * t706 + t709 * t817) * t843) + t861 * t709) * t642 + (t785 * (-(t668 * t708 - t705 * t814) * t838 + (t644 * t708 + t705 * t662) * t734 + (t723 * t705 + t708 * t817) * t844) + t862 * t708) * t641; (t783 * ((t670 * t707 + t710 * t810) * t836 + (-t646 * t707 + t710 * t664) * t738 + (-t707 * t817 + t710 * t723) * t842) - t860 * t707) * t643 + (t784 * ((t669 * t706 + t709 * t812) * t837 + (-t645 * t706 + t709 * t663) * t736 + (-t706 * t817 + t709 * t723) * t843) - t861 * t706) * t642 + (t785 * ((t668 * t705 + t708 * t814) * t838 + (-t644 * t705 + t708 * t662) * t734 + (-t705 * t817 + t708 * t723) * t844) - t862 * t705) * t641; (-t764 * t649 + t640 * t758 + t783 * (-t673 * t836 - t682 * t720 * t738 + (pkin(2) * t811 + t738 * t752) * t722)) * t643 + (-t765 * t648 + t639 * t759 + t784 * (-t672 * t837 - t681 * t720 * t736 + (pkin(2) * t813 + t736 * t753) * t722)) * t642 + (-t766 * t647 + t638 * t760 + t785 * (-t671 * t838 - t680 * t720 * t734 + (pkin(2) * t815 + t734 * t754) * t722)) * t641;];
taucX  = t1;
