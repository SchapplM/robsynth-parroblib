% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRRR8V2G2A0
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
% Datum: 2020-08-06 17:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G2A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:49:29
% EndTime: 2020-08-06 17:49:35
% DurationCPUTime: 5.49s
% Computational Cost: add. (33657->292), mult. (75984->556), div. (4716->7), fcn. (74370->22), ass. (0->208)
t732 = sin(qJ(2,1));
t738 = cos(qJ(2,1));
t744 = pkin(7) + pkin(6);
t678 = pkin(2) * t732 - t738 * t744;
t720 = sin(pkin(4));
t722 = cos(pkin(4));
t731 = sin(qJ(3,1));
t797 = t722 * t731;
t663 = pkin(3) * t797 + t678 * t720;
t737 = cos(qJ(3,1));
t812 = t720 * t732;
t718 = t737 ^ 2;
t841 = pkin(3) * t718;
t642 = 0.1e1 / (pkin(2) * t797 + t663 * t737 + t812 * t841);
t730 = sin(qJ(2,2));
t736 = cos(qJ(2,2));
t677 = pkin(2) * t730 - t736 * t744;
t729 = sin(qJ(3,2));
t799 = t722 * t729;
t662 = pkin(3) * t799 + t677 * t720;
t735 = cos(qJ(3,2));
t814 = t720 * t730;
t717 = t735 ^ 2;
t842 = pkin(3) * t717;
t641 = 0.1e1 / (pkin(2) * t799 + t662 * t735 + t814 * t842);
t728 = sin(qJ(2,3));
t734 = cos(qJ(2,3));
t676 = pkin(2) * t728 - t734 * t744;
t727 = sin(qJ(3,3));
t801 = t722 * t727;
t661 = pkin(3) * t801 + t676 * t720;
t733 = cos(qJ(3,3));
t816 = t720 * t728;
t716 = t733 ^ 2;
t843 = pkin(3) * t716;
t640 = 0.1e1 / (pkin(2) * t801 + t661 * t733 + t816 * t843);
t679 = pkin(2) * t734 + t728 * t744;
t719 = sin(pkin(8));
t721 = cos(pkin(8));
t795 = t722 * t734;
t819 = t719 * t722;
t840 = pkin(3) * t733;
t637 = (t719 * t795 + t721 * t728) * t840 + t679 * t819 + t676 * t721;
t800 = t722 * t728;
t667 = t719 * t800 - t721 * t734;
t811 = t720 * t733;
t646 = t667 * t727 + t719 * t811;
t670 = t719 * t734 + t721 * t800;
t804 = t721 * t733;
t649 = -t670 * t727 - t720 * t804;
t741 = xDP(3);
t724 = legFrame(3,2);
t704 = sin(t724);
t707 = cos(t724);
t742 = xDP(2);
t743 = xDP(1);
t756 = t704 * t742 - t707 * t743;
t622 = (t756 * t646 + t649 * t741) * t640;
t827 = t622 * t744;
t771 = t727 * t827;
t805 = t721 * t722;
t634 = (t719 * t728 - t721 * t795) * t840 - t679 * t805 + t719 * t676;
t748 = 0.1e1 / pkin(3);
t616 = (t634 * t741 + t756 * t637) * t748 * t640;
t846 = pkin(3) * t616;
t607 = t771 - t846;
t774 = t727 * t846;
t788 = t728 * t733;
t810 = t720 * t734;
t817 = t720 * t727;
t828 = t622 * t734;
t589 = (((t616 * t722 + t622 * t810) * t843 + ((-t774 + t827) * t728 + pkin(2) * t828) * t811 + t722 * t607) * t622 + (t616 * t810 + (t716 * t722 - t788 * t817 - t722) * t622) * t846) * t640;
t664 = pkin(3) * t788 + t676;
t792 = t722 * t748;
t749 = pkin(2) ^ 2;
t697 = t744 ^ 2 + t749;
t747 = pkin(3) ^ 2;
t785 = 0.2e1 * pkin(2) * pkin(3);
t831 = (-t744 * t774 + (t716 * t747 + t733 * t785 + t697) * t622) * t622;
t592 = t640 * t792 * t831 + (-t722 * t771 + (-t664 * t817 + t722 * (pkin(2) * t733 + t843)) * t616) / (t664 * t811 + (pkin(2) + t840) * t801) * t616;
t598 = (-t733 * t831 - (pkin(2) * t616 - t607 * t733) * t846) * t640;
t619 = t622 ^ 2;
t837 = mrSges(3,1) * t727;
t685 = mrSges(3,2) * t733 + t837;
t768 = t733 * mrSges(3,1) - mrSges(3,2) * t727;
t652 = t685 * t816 - t768 * t722;
t699 = mrSges(3,2) * pkin(6) - Ifges(3,6);
t700 = mrSges(3,1) * pkin(6) - Ifges(3,5);
t673 = t699 * t733 + t700 * t727;
t723 = Ifges(3,1) - Ifges(3,2);
t739 = pkin(2) * mrSges(3,2);
t852 = -2 * Ifges(3,4);
t853 = Ifges(3,4) + t716 * t852 + (-t723 * t727 + t739) * t733;
t759 = (-Ifges(3,3) * t592 + t589 * t673 + t598 * t652 + t619 * (pkin(2) * t837 + t853)) * t748;
t740 = mrSges(3,1) * pkin(2);
t750 = -0.2e1 * pkin(6) * mrSges(3,3) + (-pkin(6) ^ 2 - t749) * m(3) - Ifges(3,1) - Ifges(2,3);
t778 = 0.2e1 * t739;
t698 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t710 = m(3) * pkin(2) + mrSges(2,1);
t822 = ((t768 + t710) * t734 + t728 * t698) * t720;
t850 = t700 / 0.2e1;
t851 = -t699 / 0.2e1;
t765 = 0.2e1 * ((t851 * t727 + t850 * t733) * t616 + (t740 * t727 + t853) * t622) * t616 + t598 * t822 - (t723 * t716 - 0.2e1 * (Ifges(3,4) * t727 + t740) * t733 + t727 * t778 + t750) * t589 - t673 * t592;
t867 = -t637 * t759 + t765 * t646;
t680 = pkin(2) * t736 + t730 * t744;
t794 = t722 * t736;
t839 = pkin(3) * t735;
t638 = (t719 * t794 + t721 * t730) * t839 + t680 * t819 + t677 * t721;
t798 = t722 * t730;
t668 = t719 * t798 - t721 * t736;
t809 = t720 * t735;
t647 = t668 * t729 + t719 * t809;
t671 = t719 * t736 + t721 * t798;
t803 = t721 * t735;
t650 = -t671 * t729 - t720 * t803;
t725 = legFrame(2,2);
t705 = sin(t725);
t708 = cos(t725);
t755 = t705 * t742 - t708 * t743;
t623 = (t755 * t647 + t650 * t741) * t641;
t825 = t623 * t744;
t770 = t729 * t825;
t635 = (t719 * t730 - t721 * t794) * t839 - t680 * t805 + t719 * t677;
t617 = (t635 * t741 + t755 * t638) * t748 * t641;
t845 = pkin(3) * t617;
t608 = t770 - t845;
t773 = t729 * t845;
t787 = t730 * t735;
t808 = t720 * t736;
t815 = t720 * t729;
t826 = t623 * t736;
t590 = (((t617 * t722 + t623 * t808) * t842 + ((-t773 + t825) * t730 + pkin(2) * t826) * t809 + t722 * t608) * t623 + (t617 * t808 + (t717 * t722 - t787 * t815 - t722) * t623) * t845) * t641;
t665 = pkin(3) * t787 + t677;
t830 = (-t744 * t773 + (t717 * t747 + t735 * t785 + t697) * t623) * t623;
t593 = t641 * t792 * t830 + (-t722 * t770 + (-t665 * t815 + t722 * (pkin(2) * t735 + t842)) * t617) / (t665 * t809 + (pkin(2) + t839) * t799) * t617;
t599 = (-t735 * t830 - (pkin(2) * t617 - t608 * t735) * t845) * t641;
t620 = t623 ^ 2;
t836 = mrSges(3,1) * t729;
t686 = mrSges(3,2) * t735 + t836;
t767 = t735 * mrSges(3,1) - mrSges(3,2) * t729;
t653 = t686 * t814 - t767 * t722;
t674 = t699 * t735 + t700 * t729;
t854 = Ifges(3,4) + t717 * t852 + (-t723 * t729 + t739) * t735;
t758 = (-Ifges(3,3) * t593 + t590 * t674 + t599 * t653 + t620 * (pkin(2) * t836 + t854)) * t748;
t821 = ((t767 + t710) * t736 + t730 * t698) * t720;
t764 = 0.2e1 * ((t851 * t729 + t850 * t735) * t617 + (t740 * t729 + t854) * t623) * t617 + t599 * t821 - (t723 * t717 - 0.2e1 * (Ifges(3,4) * t729 + t740) * t735 + t729 * t778 + t750) * t590 - t674 * t593;
t866 = -t638 * t758 + t764 * t647;
t681 = pkin(2) * t738 + t732 * t744;
t793 = t722 * t738;
t838 = pkin(3) * t737;
t639 = (t719 * t793 + t721 * t732) * t838 + t681 * t819 + t678 * t721;
t796 = t722 * t732;
t669 = t719 * t796 - t721 * t738;
t807 = t720 * t737;
t648 = t669 * t731 + t719 * t807;
t672 = t719 * t738 + t721 * t796;
t802 = t721 * t737;
t651 = -t672 * t731 - t720 * t802;
t726 = legFrame(1,2);
t706 = sin(t726);
t709 = cos(t726);
t754 = t706 * t742 - t709 * t743;
t624 = (t754 * t648 + t651 * t741) * t642;
t823 = t624 * t744;
t769 = t731 * t823;
t636 = (t719 * t732 - t721 * t793) * t838 - t681 * t805 + t719 * t678;
t618 = (t636 * t741 + t754 * t639) * t748 * t642;
t844 = pkin(3) * t618;
t609 = t769 - t844;
t772 = t731 * t844;
t786 = t732 * t737;
t806 = t720 * t738;
t813 = t720 * t731;
t824 = t624 * t738;
t591 = (((t618 * t722 + t624 * t806) * t841 + ((-t772 + t823) * t732 + pkin(2) * t824) * t807 + t722 * t609) * t624 + (t618 * t806 + (t718 * t722 - t786 * t813 - t722) * t624) * t844) * t642;
t666 = pkin(3) * t786 + t678;
t829 = (-t744 * t772 + (t718 * t747 + t737 * t785 + t697) * t624) * t624;
t594 = t642 * t792 * t829 + (-t722 * t769 + (-t666 * t813 + t722 * (pkin(2) * t737 + t841)) * t618) / (t666 * t807 + (pkin(2) + t838) * t797) * t618;
t600 = (-t737 * t829 - (pkin(2) * t618 - t609 * t737) * t844) * t642;
t621 = t624 ^ 2;
t835 = mrSges(3,1) * t731;
t687 = mrSges(3,2) * t737 + t835;
t766 = t737 * mrSges(3,1) - mrSges(3,2) * t731;
t654 = t687 * t812 - t766 * t722;
t675 = t699 * t737 + t700 * t731;
t855 = Ifges(3,4) + t718 * t852 + (-t723 * t731 + t739) * t737;
t757 = (-Ifges(3,3) * t594 + t591 * t675 + t600 * t654 + t621 * (pkin(2) * t835 + t855)) * t748;
t820 = ((t766 + t710) * t738 + t732 * t698) * t720;
t763 = 0.2e1 * ((t851 * t731 + t850 * t737) * t618 + (t740 * t731 + t855) * t624) * t618 + t600 * t820 - (t723 * t718 - 0.2e1 * (Ifges(3,4) * t731 + t740) * t737 + t731 * t778 + t750) * t591 - t675 * t594;
t865 = -t639 * t757 + t763 * t648;
t849 = pkin(2) * t727;
t848 = pkin(2) * t729;
t847 = pkin(2) * t731;
t818 = t720 * t721;
t613 = t616 ^ 2;
t714 = -m(1) - m(2) - m(3);
t784 = -t589 * t822 + t592 * t652 + t598 * t714 + ((-t619 * t710 - t768 * (t619 + t613)) * t728 + (-0.2e1 * t616 * t685 + t622 * t698) * t828) * t720 - t613 * t722 * t685;
t614 = t617 ^ 2;
t783 = -t590 * t821 + t593 * t653 + t599 * t714 + ((-t620 * t710 - t767 * (t620 + t614)) * t730 + (-0.2e1 * t617 * t686 + t623 * t698) * t826) * t720 - t614 * t722 * t686;
t615 = t618 ^ 2;
t782 = -t591 * t820 + t594 * t654 + t600 * t714 + ((-t621 * t710 - t766 * (t621 + t615)) * t732 + (-0.2e1 * t618 * t687 + t624 * t698) * t824) * t720 - t615 * t722 * t687;
t753 = pkin(3) * t817 - t676 * t722;
t752 = pkin(3) * t815 - t677 * t722;
t751 = pkin(3) * t813 - t678 * t722;
t645 = -t681 * t719 + t751 * t721;
t644 = -t680 * t719 + t752 * t721;
t643 = -t679 * t719 + t753 * t721;
t1 = [(t782 * ((t672 * t709 + t706 * t812) * t841 + (-t645 * t709 + t663 * t706) * t737 + (t706 * t722 - t709 * t818) * t847) + t865 * t709) * t642 + (t783 * ((t671 * t708 + t705 * t814) * t842 + (-t644 * t708 + t662 * t705) * t735 + (t705 * t722 - t708 * t818) * t848) + t866 * t708) * t641 + (t784 * ((t670 * t707 + t704 * t816) * t843 + (-t643 * t707 + t661 * t704) * t733 + (t704 * t722 - t707 * t818) * t849) + t867 * t707) * t640; (t782 * (-(t672 * t706 - t709 * t812) * t841 + (t645 * t706 + t663 * t709) * t737 + (t706 * t818 + t709 * t722) * t847) - t865 * t706) * t642 + (t783 * (-(t671 * t705 - t708 * t814) * t842 + (t644 * t705 + t662 * t708) * t735 + (t705 * t818 + t708 * t722) * t848) - t866 * t705) * t641 + (t784 * (-(t670 * t704 - t707 * t816) * t843 + (t643 * t704 + t661 * t707) * t733 + (t704 * t818 + t707 * t722) * t849) - t867 * t704) * t640; (-t763 * t651 + t636 * t757 + t782 * (-t669 * t841 + t681 * t802 + (pkin(2) * t813 + t751 * t737) * t719)) * t642 + (-t764 * t650 + t635 * t758 + t783 * (-t668 * t842 + t680 * t803 + (pkin(2) * t815 + t752 * t735) * t719)) * t641 + (-t765 * t649 + t634 * t759 + t784 * (-t667 * t843 + t679 * t804 + (pkin(2) * t817 + t753 * t733) * t719)) * t640;];
taucX  = t1;
