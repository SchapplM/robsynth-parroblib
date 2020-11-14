% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRRR8V1G3A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
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
% Datum: 2020-08-06 17:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:15:56
% EndTime: 2020-08-06 17:16:00
% DurationCPUTime: 3.61s
% Computational Cost: add. (15369->231), mult. (49278->489), div. (6885->8), fcn. (54219->28), ass. (0->208)
t747 = cos(qJ(3,2));
t748 = cos(qJ(2,2));
t742 = sin(qJ(2,2));
t812 = t742 * t747;
t701 = pkin(2) * t812 - t748 * pkin(5);
t733 = sin(pkin(3));
t735 = cos(pkin(3));
t741 = sin(qJ(3,2));
t845 = pkin(2) * t741;
t695 = t701 * t733 + t735 * t845;
t862 = 0.1e1 / t695;
t831 = t862 / t747;
t745 = cos(qJ(3,3));
t746 = cos(qJ(2,3));
t740 = sin(qJ(2,3));
t814 = t740 * t745;
t700 = pkin(2) * t814 - t746 * pkin(5);
t739 = sin(qJ(3,3));
t846 = pkin(2) * t739;
t694 = t700 * t733 + t735 * t846;
t863 = 0.1e1 / t694;
t832 = t863 / t745;
t743 = sin(qJ(3,1));
t749 = cos(qJ(3,1));
t781 = rSges(3,1) * t749 - rSges(3,2) * t743;
t867 = t781 * m(3);
t782 = rSges(3,1) * t747 - rSges(3,2) * t741;
t866 = t782 * m(3);
t783 = rSges(3,1) * t745 - rSges(3,2) * t739;
t865 = t783 * m(3);
t744 = sin(qJ(2,1));
t810 = t744 * t749;
t793 = t733 * t810;
t821 = t735 * t743;
t750 = cos(qJ(2,1));
t824 = t733 * t750;
t864 = 0.1e1 / (-pkin(5) * t824 + (t793 + t821) * pkin(2));
t718 = m(3) * rSges(3,1) * rSges(3,2) - Icges(3,4);
t861 = 0.2e1 * t718;
t753 = m(2) * rSges(2,1);
t860 = m(3) * rSges(3,3);
t760 = rSges(3,2) ^ 2;
t761 = rSges(3,1) ^ 2;
t710 = (-t760 + t761) * m(3) - Icges(3,1) + Icges(3,2);
t859 = -t710 / 0.2e1;
t716 = rSges(3,2) * t860 - Icges(3,6);
t858 = -t716 / 0.4e1;
t717 = rSges(3,1) * t860 - Icges(3,5);
t857 = t717 / 0.4e1;
t706 = t739 * rSges(3,1) + t745 * rSges(3,2);
t856 = m(3) * (-t733 * t740 * t706 + t783 * t735);
t707 = t741 * rSges(3,1) + t747 * rSges(3,2);
t855 = m(3) * (-t733 * t742 * t707 + t782 * t735);
t708 = t743 * rSges(3,1) + t749 * rSges(3,2);
t854 = m(3) * (-t733 * t744 * t708 + t781 * t735);
t853 = m(3) * t735;
t732 = sin(pkin(6));
t734 = cos(pkin(6));
t819 = t735 * t746;
t823 = t735 * t740;
t844 = pkin(2) * t745;
t665 = (-t732 * t740 + t734 * t819) * t844 + pkin(5) * (t732 * t746 + t734 * t823);
t668 = (t732 * t819 + t734 * t740) * t844 + (t732 * t823 - t734 * t746) * pkin(5);
t754 = xDP(3);
t764 = 0.1e1 / pkin(2);
t736 = legFrame(3,2);
t719 = sin(t736);
t722 = cos(t736);
t755 = xDP(2);
t756 = xDP(1);
t779 = t719 * t755 - t722 * t756;
t644 = (t779 * t665 + t668 * t754) * t764 * t832;
t852 = pkin(2) * t644;
t818 = t735 * t748;
t822 = t735 * t742;
t843 = pkin(2) * t747;
t666 = (-t732 * t742 + t734 * t818) * t843 + pkin(5) * (t732 * t748 + t734 * t822);
t669 = (t732 * t818 + t734 * t742) * t843 + (t732 * t822 - t734 * t748) * pkin(5);
t737 = legFrame(2,2);
t720 = sin(t737);
t723 = cos(t737);
t778 = t720 * t755 - t723 * t756;
t645 = (t778 * t666 + t669 * t754) * t764 * t831;
t851 = pkin(2) * t645;
t817 = t735 * t750;
t820 = t735 * t744;
t842 = pkin(2) * t749;
t667 = (-t732 * t744 + t734 * t817) * t842 + pkin(5) * (t732 * t750 + t734 * t820);
t670 = (t732 * t817 + t734 * t744) * t842 + (t732 * t820 - t734 * t750) * pkin(5);
t738 = legFrame(1,2);
t721 = sin(t738);
t724 = cos(t738);
t777 = t721 * t755 - t724 * t756;
t731 = 0.1e1 / t749;
t830 = t864 * t731;
t646 = (t777 * t667 + t670 * t754) * t764 * t830;
t850 = pkin(2) * t646;
t726 = t745 ^ 2;
t849 = pkin(2) * t726;
t728 = t747 ^ 2;
t848 = pkin(2) * t728;
t730 = t749 ^ 2;
t847 = pkin(2) * t730;
t767 = t733 * t745 + t739 * t823;
t815 = t739 * t746;
t674 = t767 * t732 - t734 * t815;
t677 = t732 * t815 + t767 * t734;
t650 = (t674 * t754 + t779 * t677) * t832;
t841 = pkin(5) * t650;
t766 = t733 * t747 + t741 * t822;
t813 = t741 * t748;
t675 = t766 * t732 - t734 * t813;
t678 = t732 * t813 + t766 * t734;
t651 = (t675 * t754 + t778 * t678) * t831;
t840 = pkin(5) * t651;
t765 = t733 * t749 + t743 * t820;
t811 = t743 * t750;
t676 = t765 * t732 - t734 * t811;
t679 = t732 * t811 + t765 * t734;
t652 = (t676 * t754 + t777 * t679) * t830;
t839 = pkin(5) * t652;
t762 = pkin(5) ^ 2;
t763 = pkin(2) ^ 2;
t807 = t644 * t846;
t838 = (-pkin(5) * t807 + (t726 * t763 + t762) * t650) * t650;
t806 = t645 * t845;
t837 = (-pkin(5) * t806 + (t728 * t763 + t762) * t651) * t651;
t836 = t652 * t864;
t715 = m(2) * rSges(2,2) - t860;
t835 = ((t753 + t865) * t746 - t715 * t740) * t733;
t834 = ((t753 + t866) * t748 - t715 * t742) * t733;
t833 = ((t753 + t867) * t750 - t715 * t744) * t733;
t829 = t733 * t739;
t828 = t733 * t741;
t827 = t733 * t743;
t826 = t733 * t746;
t825 = t733 * t748;
t816 = t735 * t764;
t809 = t760 + t761;
t808 = 0.2e1 * m(3);
t805 = t743 * t850;
t804 = t739 * t841;
t803 = t741 * t840;
t802 = t743 * t839;
t801 = t722 * t832;
t800 = t723 * t831;
t799 = t724 * t830;
t798 = t710 * t739 * t745;
t797 = t710 * t741 * t747;
t796 = t710 * t743 * t749;
t795 = t733 * t814;
t794 = t733 * t812;
t635 = t804 - t852;
t617 = (((t735 * t644 + t650 * t826) * t849 - (t807 - t841) * t795 + t735 * t635) * t650 - (-t644 * t826 + (-t726 * t735 + t739 * t795 + t735) * t650) * t852) * t832;
t623 = (t816 * t838 + (-t644 * t700 * t829 + t735 * (t644 * t849 - t804)) * t644) * t832;
t626 = (t635 * t852 - t838) * t863;
t641 = t644 ^ 2;
t647 = t650 ^ 2;
t725 = -m(1) - m(2) - m(3);
t792 = (-t617 * t835 - t623 * t856 + t725 * t626 + ((-t647 * t753 - (t647 + t641) * t865) * t740 - t746 * (t706 * t644 * t808 + t650 * t715) * t650) * t733 - t641 * t706 * t853) * t863;
t636 = t803 - t851;
t618 = (((t735 * t645 + t651 * t825) * t848 - (t806 - t840) * t794 + t735 * t636) * t651 + (t645 * t825 + (t728 * t735 - t741 * t794 - t735) * t651) * t851) * t831;
t624 = (t816 * t837 + (-t645 * t701 * t828 + t735 * (t645 * t848 - t803)) * t645) * t831;
t627 = (t636 * t851 - t837) * t862;
t642 = t645 ^ 2;
t648 = t651 ^ 2;
t791 = (-t618 * t834 - t624 * t855 + t725 * t627 + ((-t648 * t753 - (t648 + t642) * t866) * t742 - t748 * (t707 * t645 * t808 + t651 * t715) * t651) * t733 - t642 * t707 * t853) * t862;
t637 = t802 - t850;
t619 = (((t735 * t646 + t652 * t824) * t847 - (t805 - t839) * t793 + t735 * t637) * t836 + (t646 * t824 + (t730 * t735 - t743 * t793 - t735) * t652) * t864 * t850) * t731;
t634 = -pkin(5) * t805 + (t730 * t763 + t762) * t652;
t702 = pkin(2) * t810 - t750 * pkin(5);
t696 = pkin(2) * t821 + t702 * t733;
t693 = 0.1e1 / t696;
t625 = (t634 * t816 * t836 + (-t646 * t702 * t827 + t735 * (t646 * t847 - t802)) * t693 * t646) * t731;
t628 = (t634 * t652 - t637 * t850) * t864;
t643 = t646 ^ 2;
t649 = t652 ^ 2;
t790 = (-t619 * t833 - t625 * t854 - t725 * t628 + ((-t649 * t753 - (t649 + t643) * t867) * t744 - t750 * (t708 * t646 * t808 + t652 * t715) * t652) * t733 - t643 * t708 * t853) * t693;
t697 = t716 * t745 + t717 * t739;
t757 = 0.2e1 * qJ(3,3);
t780 = -(rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) - (0.2e1 * rSges(3,3) ^ 2 + t809) * m(3) / 0.2e1 - Icges(3,2) / 0.2e1 - Icges(3,1) / 0.2e1 - Icges(2,3);
t789 = 0.4e1 * ((t739 * t858 + t745 * t857) * t644 + (t798 / 0.2e1 + (t726 - 0.1e1 / 0.2e1) * t718) * t650) * t644 + t626 * t835 - (cos(t757) * t859 + t718 * sin(t757) + t780) * t617 - t697 * t623;
t698 = t716 * t747 + t717 * t741;
t758 = 0.2e1 * qJ(3,2);
t788 = 0.4e1 * ((t741 * t858 + t747 * t857) * t645 + (t797 / 0.2e1 + (t728 - 0.1e1 / 0.2e1) * t718) * t651) * t645 + t627 * t834 - (cos(t758) * t859 + t718 * sin(t758) + t780) * t618 - t698 * t624;
t699 = t716 * t749 + t717 * t743;
t759 = 0.2e1 * qJ(3,1);
t787 = 0.4e1 * ((t743 * t858 + t749 * t857) * t646 + (t796 / 0.2e1 + (t730 - 0.1e1 / 0.2e1) * t718) * t652) * t646 - t628 * t833 - (cos(t759) * t859 + t718 * sin(t759) + t780) * t619 - t699 * t625;
t712 = -t809 * m(3) - Icges(3,3);
t786 = t697 * t617 + t712 * t623 - t626 * t856 + t647 * (t726 * t861 - t718 + t798);
t785 = t698 * t618 + t712 * t624 - t627 * t855 + t648 * (t728 * t861 - t718 + t797);
t784 = t699 * t619 + t712 * t625 + t628 * t854 + t649 * (t730 * t861 - t718 + t796);
t776 = pkin(2) * t829 - t700 * t735;
t775 = pkin(2) * t828 - t701 * t735;
t774 = pkin(2) * t827 - t702 * t735;
t773 = t789 * t832;
t772 = t786 * t832;
t771 = t788 * t831;
t770 = t785 * t831;
t769 = t787 * t830;
t768 = t784 * t830;
t705 = pkin(5) * t744 + t750 * t842;
t704 = pkin(5) * t742 + t748 * t843;
t703 = pkin(5) * t740 + t746 * t844;
t661 = t734 * t705 + t774 * t732;
t660 = t734 * t704 + t775 * t732;
t659 = t734 * t703 + t776 * t732;
t1 = [t787 * t679 * t799 + t788 * t678 * t800 + t789 * t677 * t801 + (t661 * t724 + t721 * t696) * t790 + (t660 * t723 + t720 * t695) * t791 + (t659 * t722 + t719 * t694) * t792 + (-t786 * t665 * t801 - t785 * t666 * t800 - t784 * t667 * t799) * t764; -t721 * t679 * t769 - t720 * t678 * t771 - t719 * t677 * t773 + (-t661 * t721 + t724 * t696) * t790 + (-t660 * t720 + t723 * t695) * t791 + (-t659 * t719 + t722 * t694) * t792 + (t719 * t665 * t772 + t720 * t666 * t770 + t721 * t667 * t768) * t764; -t676 * t769 - t675 * t771 - t674 * t773 + (-t732 * t705 + t774 * t734) * t790 + (-t732 * t704 + t775 * t734) * t791 + (-t732 * t703 + t776 * t734) * t792 + (t668 * t772 + t669 * t770 + t670 * t768) * t764;];
taucX  = t1;
