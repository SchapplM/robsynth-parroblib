% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRRR8V1G2A0
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
% Datum: 2020-08-06 17:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRRR8V1G2A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G2A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V1G2A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G2A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G2A0_coriolisvec_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G2A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V1G2A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR8V1G2A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G2A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G2A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:02:37
% EndTime: 2020-08-06 17:02:40
% DurationCPUTime: 3.62s
% Computational Cost: add. (15369->233), mult. (49278->496), div. (6885->9), fcn. (54219->28), ass. (0->211)
t749 = cos(qJ(3,3));
t750 = cos(qJ(2,3));
t744 = sin(qJ(2,3));
t818 = t744 * t749;
t705 = pkin(2) * t818 - t750 * pkin(5);
t737 = sin(pkin(3));
t739 = cos(pkin(3));
t743 = sin(qJ(3,3));
t850 = pkin(2) * t743;
t699 = t705 * t737 + t739 * t850;
t866 = 0.1e1 / t699;
t837 = t866 / t749;
t747 = sin(qJ(3,1));
t753 = cos(qJ(3,1));
t785 = rSges(3,1) * t753 - rSges(3,2) * t747;
t871 = t785 * m(3);
t745 = sin(qJ(3,2));
t751 = cos(qJ(3,2));
t786 = rSges(3,1) * t751 - rSges(3,2) * t745;
t870 = t786 * m(3);
t787 = rSges(3,1) * t749 - rSges(3,2) * t743;
t869 = t787 * m(3);
t746 = sin(qJ(2,2));
t816 = t746 * t751;
t801 = t737 * t816;
t827 = t739 * t745;
t752 = cos(qJ(2,2));
t830 = t737 * t752;
t868 = 0.1e1 / (-pkin(5) * t830 + (t801 + t827) * pkin(2));
t748 = sin(qJ(2,1));
t814 = t748 * t753;
t800 = t737 * t814;
t825 = t739 * t747;
t754 = cos(qJ(2,1));
t829 = t737 * t754;
t867 = 0.1e1 / (-pkin(5) * t829 + (t800 + t825) * pkin(2));
t722 = m(3) * rSges(3,1) * rSges(3,2) - Icges(3,4);
t865 = 0.2e1 * t722;
t757 = m(2) * rSges(2,1);
t864 = m(3) * rSges(3,3);
t764 = rSges(3,2) ^ 2;
t765 = rSges(3,1) ^ 2;
t715 = (-t764 + t765) * m(3) - Icges(3,1) + Icges(3,2);
t863 = -t715 / 0.2e1;
t720 = rSges(3,2) * t864 - Icges(3,6);
t862 = -t720 / 0.4e1;
t721 = rSges(3,1) * t864 - Icges(3,5);
t861 = t721 / 0.4e1;
t711 = t743 * rSges(3,1) + t749 * rSges(3,2);
t860 = m(3) * (-t737 * t744 * t711 + t787 * t739);
t712 = t745 * rSges(3,1) + t751 * rSges(3,2);
t859 = m(3) * (-t737 * t746 * t712 + t786 * t739);
t713 = t747 * rSges(3,1) + t753 * rSges(3,2);
t858 = m(3) * (-t737 * t748 * t713 + t785 * t739);
t857 = m(3) * t739;
t736 = sin(pkin(6));
t738 = cos(pkin(6));
t823 = t739 * t750;
t828 = t739 * t744;
t849 = pkin(2) * t749;
t668 = -(-t736 * t744 + t738 * t823) * t849 - pkin(5) * (t736 * t750 + t738 * t828);
t671 = (t736 * t823 + t738 * t744) * t849 + (t736 * t828 - t738 * t750) * pkin(5);
t758 = xDP(3);
t768 = 0.1e1 / pkin(2);
t740 = legFrame(3,2);
t723 = sin(t740);
t726 = cos(t740);
t759 = xDP(2);
t760 = xDP(1);
t783 = t723 * t759 - t726 * t760;
t647 = (t668 * t758 + t783 * t671) * t768 * t837;
t856 = pkin(2) * t647;
t822 = t739 * t752;
t826 = t739 * t746;
t848 = pkin(2) * t751;
t669 = -(-t736 * t746 + t738 * t822) * t848 - pkin(5) * (t736 * t752 + t738 * t826);
t672 = (t736 * t822 + t738 * t746) * t848 + (t736 * t826 - t738 * t752) * pkin(5);
t741 = legFrame(2,2);
t724 = sin(t741);
t727 = cos(t741);
t782 = t724 * t759 - t727 * t760;
t733 = 0.1e1 / t751;
t836 = t868 * t733;
t648 = (t669 * t758 + t782 * t672) * t768 * t836;
t855 = pkin(2) * t648;
t821 = t739 * t754;
t824 = t739 * t748;
t847 = pkin(2) * t753;
t670 = -(-t736 * t748 + t738 * t821) * t847 - pkin(5) * (t736 * t754 + t738 * t824);
t673 = (t736 * t821 + t738 * t748) * t847 + (t736 * t824 - t738 * t754) * pkin(5);
t742 = legFrame(1,2);
t725 = sin(t742);
t728 = cos(t742);
t781 = t725 * t759 - t728 * t760;
t735 = 0.1e1 / t753;
t835 = t867 * t735;
t649 = (t670 * t758 + t781 * t673) * t768 * t835;
t854 = pkin(2) * t649;
t730 = t749 ^ 2;
t853 = pkin(2) * t730;
t732 = t751 ^ 2;
t852 = pkin(2) * t732;
t734 = t753 ^ 2;
t851 = pkin(2) * t734;
t771 = t737 * t749 + t743 * t828;
t819 = t743 * t750;
t677 = t771 * t736 - t738 * t819;
t680 = -t736 * t819 - t771 * t738;
t653 = (t783 * t677 + t680 * t758) * t837;
t846 = pkin(5) * t653;
t770 = t737 * t751 + t745 * t826;
t817 = t745 * t752;
t678 = t770 * t736 - t738 * t817;
t681 = -t736 * t817 - t770 * t738;
t654 = (t782 * t678 + t681 * t758) * t836;
t845 = pkin(5) * t654;
t769 = t737 * t753 + t747 * t824;
t815 = t747 * t754;
t679 = t769 * t736 - t738 * t815;
t682 = -t736 * t815 - t769 * t738;
t655 = (t781 * t679 + t682 * t758) * t835;
t844 = pkin(5) * t655;
t766 = pkin(5) ^ 2;
t767 = pkin(2) ^ 2;
t811 = t647 * t850;
t843 = (-pkin(5) * t811 + (t730 * t767 + t766) * t653) * t653;
t842 = t654 * t868;
t841 = t655 * t867;
t719 = m(2) * rSges(2,2) - t864;
t840 = ((t757 + t869) * t750 - t719 * t744) * t737;
t839 = ((t757 + t870) * t752 - t719 * t746) * t737;
t838 = ((t757 + t871) * t754 - t719 * t748) * t737;
t834 = t737 * t743;
t833 = t737 * t745;
t832 = t737 * t747;
t831 = t737 * t750;
t820 = t739 * t768;
t813 = t764 + t765;
t812 = 0.2e1 * m(3);
t810 = t745 * t855;
t809 = t747 * t854;
t808 = t743 * t846;
t807 = t745 * t845;
t806 = t747 * t844;
t805 = t726 * t837;
t804 = t727 * t836;
t803 = t728 * t835;
t802 = t737 * t818;
t799 = t743 * t715 * t749;
t798 = t745 * t715 * t751;
t797 = t747 * t715 * t753;
t638 = t808 - t856;
t620 = (((t739 * t647 + t653 * t831) * t853 - (t811 - t846) * t802 + t739 * t638) * t653 - (-t647 * t831 + (-t730 * t739 + t743 * t802 + t739) * t653) * t856) * t837;
t626 = (t820 * t843 + (-t647 * t705 * t834 + t739 * (t647 * t853 - t808)) * t647) * t837;
t629 = (t638 * t856 - t843) * t866;
t644 = t647 ^ 2;
t650 = t653 ^ 2;
t729 = -m(1) - m(2) - m(3);
t796 = (-t620 * t840 - t626 * t860 + t729 * t629 + ((-t650 * t757 - (t650 + t644) * t869) * t744 - (t711 * t647 * t812 + t653 * t719) * t653 * t750) * t737 - t644 * t711 * t857) * t866;
t639 = t807 - t855;
t621 = (((t739 * t648 + t654 * t830) * t852 - (t810 - t845) * t801 + t739 * t639) * t842 + (t648 * t830 + (t732 * t739 - t745 * t801 - t739) * t654) * t868 * t855) * t733;
t636 = -pkin(5) * t810 + (t732 * t767 + t766) * t654;
t706 = pkin(2) * t816 - t752 * pkin(5);
t700 = pkin(2) * t827 + t706 * t737;
t697 = 0.1e1 / t700;
t627 = (t636 * t820 * t842 + (-t648 * t706 * t833 + t739 * (t648 * t852 - t807)) * t697 * t648) * t733;
t630 = (t636 * t654 - t639 * t855) * t868;
t645 = t648 ^ 2;
t651 = t654 ^ 2;
t795 = (-t621 * t839 - t627 * t859 - t729 * t630 + ((-t651 * t757 - (t651 + t645) * t870) * t746 - (t712 * t648 * t812 + t654 * t719) * t654 * t752) * t737 - t645 * t712 * t857) * t697;
t640 = t806 - t854;
t622 = (((t739 * t649 + t655 * t829) * t851 - (t809 - t844) * t800 + t739 * t640) * t841 + (t649 * t829 + (t734 * t739 - t747 * t800 - t739) * t655) * t867 * t854) * t735;
t637 = -pkin(5) * t809 + (t734 * t767 + t766) * t655;
t707 = pkin(2) * t814 - t754 * pkin(5);
t701 = pkin(2) * t825 + t707 * t737;
t698 = 0.1e1 / t701;
t628 = (t637 * t820 * t841 + (-t649 * t707 * t832 + t739 * (t649 * t851 - t806)) * t698 * t649) * t735;
t631 = (t637 * t655 - t640 * t854) * t867;
t646 = t649 ^ 2;
t652 = t655 ^ 2;
t794 = (-t622 * t838 - t628 * t858 - t729 * t631 + ((-t652 * t757 - (t652 + t646) * t871) * t748 - (t713 * t649 * t812 + t655 * t719) * t655 * t754) * t737 - t646 * t713 * t857) * t698;
t702 = t720 * t749 + t721 * t743;
t761 = 0.2e1 * qJ(3,3);
t784 = -(rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) - (0.2e1 * rSges(3,3) ^ 2 + t813) * m(3) / 0.2e1 - Icges(3,2) / 0.2e1 - Icges(3,1) / 0.2e1 - Icges(2,3);
t793 = 0.4e1 * ((t743 * t862 + t749 * t861) * t647 + (t799 / 0.2e1 + (t730 - 0.1e1 / 0.2e1) * t722) * t653) * t647 + t629 * t840 - (cos(t761) * t863 + t722 * sin(t761) + t784) * t620 - t702 * t626;
t703 = t720 * t751 + t721 * t745;
t762 = 0.2e1 * qJ(3,2);
t792 = 0.4e1 * ((t745 * t862 + t751 * t861) * t648 + (t798 / 0.2e1 + (t732 - 0.1e1 / 0.2e1) * t722) * t654) * t648 - t630 * t839 - (cos(t762) * t863 + t722 * sin(t762) + t784) * t621 - t703 * t627;
t704 = t720 * t753 + t721 * t747;
t763 = 0.2e1 * qJ(3,1);
t791 = 0.4e1 * ((t747 * t862 + t753 * t861) * t649 + (t797 / 0.2e1 + (t734 - 0.1e1 / 0.2e1) * t722) * t655) * t649 - t631 * t838 - (cos(t763) * t863 + t722 * sin(t763) + t784) * t622 - t704 * t628;
t717 = -t813 * m(3) - Icges(3,3);
t790 = t702 * t620 + t717 * t626 - t629 * t860 + t650 * (t730 * t865 - t722 + t799);
t789 = t703 * t621 + t717 * t627 + t630 * t859 + t651 * (t732 * t865 - t722 + t798);
t788 = t704 * t622 + t717 * t628 + t631 * t858 + t652 * (t734 * t865 - t722 + t797);
t780 = pkin(2) * t834 - t705 * t739;
t779 = pkin(2) * t833 - t706 * t739;
t778 = pkin(2) * t832 - t707 * t739;
t777 = t793 * t837;
t776 = t790 * t837;
t775 = t792 * t836;
t774 = t789 * t836;
t773 = t791 * t835;
t772 = t788 * t835;
t710 = pkin(5) * t748 + t754 * t847;
t709 = pkin(5) * t746 + t752 * t848;
t708 = pkin(5) * t744 + t750 * t849;
t664 = -t736 * t710 + t778 * t738;
t663 = -t736 * t709 + t779 * t738;
t662 = -t736 * t708 + t780 * t738;
t1 = [t791 * t679 * t803 + t792 * t678 * t804 + t793 * t677 * t805 + (-t664 * t728 + t701 * t725) * t794 + (-t663 * t727 + t700 * t724) * t795 + (-t662 * t726 + t699 * t723) * t796 + (-t790 * t671 * t805 - t789 * t672 * t804 - t788 * t673 * t803) * t768; -t725 * t679 * t773 - t724 * t678 * t775 - t723 * t677 * t777 + (t664 * t725 + t701 * t728) * t794 + (t663 * t724 + t700 * t727) * t795 + (t662 * t723 + t699 * t726) * t796 + (t723 * t671 * t776 + t724 * t672 * t774 + t725 * t673 * t772) * t768; -t682 * t773 - t681 * t775 - t680 * t777 + (t738 * t710 + t778 * t736) * t794 + (t738 * t709 + t779 * t736) * t795 + (t738 * t708 + t780 * t736) * t796 + (t668 * t776 + t669 * t774 + t670 * t772) * t768;];
taucX  = t1;
