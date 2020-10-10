% Calculate inertia matrix for parallel robot
% P4PRRRR8V2G3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [4x3]
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
% MX [4x4]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:30
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P4PRRRR8V2G3A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G3A0_inertia_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G3A0_inertia_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G3A0_inertia_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V2G3A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR8V2G3A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P4PRRRR8V2G3A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G3A0_inertia_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G3A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:25:51
% EndTime: 2020-08-07 11:25:55
% DurationCPUTime: 4.04s
% Computational Cost: add. (15484->458), mult. (32520->859), div. (1280->5), fcn. (30920->30), ass. (0->303)
t738 = sin(qJ(2,4));
t740 = cos(qJ(2,4));
t759 = pkin(7) + pkin(6);
t701 = pkin(2) * t738 - t759 * t740;
t734 = sin(pkin(4));
t736 = cos(pkin(4));
t737 = sin(qJ(3,4));
t814 = t737 * t736;
t677 = pkin(3) * t814 + t701 * t734;
t739 = cos(qJ(3,4));
t834 = t734 * t738;
t727 = t739 ^ 2;
t879 = pkin(3) * t727;
t649 = 0.1e1 / (pkin(2) * t814 + t677 * t739 + t834 * t879);
t763 = xP(4);
t725 = sin(t763);
t726 = cos(t763);
t766 = koppelP(4,2);
t770 = koppelP(4,1);
t693 = -t725 * t770 - t726 * t766;
t697 = -t725 * t766 + t726 * t770;
t742 = legFrame(4,2);
t717 = sin(t742);
t721 = cos(t742);
t887 = t649 * (t693 * t721 - t697 * t717);
t747 = sin(qJ(2,3));
t753 = cos(qJ(2,3));
t703 = pkin(2) * t747 - t759 * t753;
t746 = sin(qJ(3,3));
t813 = t746 * t736;
t678 = pkin(3) * t813 + t703 * t734;
t752 = cos(qJ(3,3));
t831 = t734 * t747;
t730 = t752 ^ 2;
t878 = pkin(3) * t730;
t650 = 0.1e1 / (pkin(2) * t813 + t678 * t752 + t831 * t878);
t767 = koppelP(3,2);
t771 = koppelP(3,1);
t694 = -t725 * t771 - t726 * t767;
t698 = -t725 * t767 + t726 * t771;
t743 = legFrame(3,2);
t718 = sin(t743);
t722 = cos(t743);
t886 = t650 * (t694 * t722 - t698 * t718);
t749 = sin(qJ(2,2));
t755 = cos(qJ(2,2));
t704 = pkin(2) * t749 - t759 * t755;
t748 = sin(qJ(3,2));
t812 = t748 * t736;
t679 = pkin(3) * t812 + t704 * t734;
t754 = cos(qJ(3,2));
t829 = t734 * t749;
t731 = t754 ^ 2;
t877 = pkin(3) * t731;
t651 = 0.1e1 / (pkin(2) * t812 + t679 * t754 + t829 * t877);
t768 = koppelP(2,2);
t772 = koppelP(2,1);
t695 = -t725 * t772 - t726 * t768;
t699 = -t725 * t768 + t726 * t772;
t744 = legFrame(2,2);
t719 = sin(t744);
t723 = cos(t744);
t885 = t651 * (t695 * t723 - t699 * t719);
t751 = sin(qJ(2,1));
t757 = cos(qJ(2,1));
t705 = pkin(2) * t751 - t759 * t757;
t750 = sin(qJ(3,1));
t811 = t750 * t736;
t680 = pkin(3) * t811 + t705 * t734;
t756 = cos(qJ(3,1));
t827 = t734 * t751;
t732 = t756 ^ 2;
t876 = pkin(3) * t732;
t652 = 0.1e1 / (pkin(2) * t811 + t680 * t756 + t827 * t876);
t769 = koppelP(1,2);
t773 = koppelP(1,1);
t696 = -t725 * t773 - t726 * t769;
t700 = -t725 * t769 + t726 * t773;
t745 = legFrame(1,2);
t720 = sin(t745);
t724 = cos(t745);
t884 = t652 * (t696 * t724 - t700 * t720);
t883 = pkin(2) * t737;
t882 = pkin(2) * t746;
t881 = pkin(2) * t748;
t880 = pkin(2) * t750;
t875 = pkin(3) * t739;
t874 = pkin(3) * t752;
t873 = pkin(3) * t754;
t872 = pkin(3) * t756;
t774 = 0.1e1 / pkin(3);
t871 = Ifges(3,3) * t774;
t870 = m(3) * pkin(2) + mrSges(2,1);
t702 = pkin(2) * t740 + t738 * t759;
t733 = sin(pkin(8));
t735 = cos(pkin(8));
t821 = t736 * t740;
t823 = t735 * t736;
t641 = (t733 * t738 - t735 * t821) * t875 - t702 * t823 + t701 * t733;
t869 = t641 * t649;
t868 = t641 * t774;
t836 = t733 * t736;
t642 = (t733 * t821 + t735 * t738) * t875 + t702 * t836 + t735 * t701;
t867 = t642 * t649;
t706 = pkin(2) * t753 + t747 * t759;
t817 = t736 * t753;
t643 = (t733 * t747 - t735 * t817) * t874 - t706 * t823 + t703 * t733;
t866 = t643 * t650;
t865 = t643 * t774;
t707 = pkin(2) * t755 + t749 * t759;
t816 = t736 * t755;
t644 = (t733 * t749 - t735 * t816) * t873 - t707 * t823 + t704 * t733;
t864 = t644 * t651;
t863 = t644 * t774;
t708 = pkin(2) * t757 + t751 * t759;
t815 = t736 * t757;
t645 = (t733 * t751 - t735 * t815) * t872 - t708 * t823 + t705 * t733;
t862 = t645 * t652;
t861 = t645 * t774;
t646 = (t733 * t817 + t735 * t747) * t874 + t706 * t836 + t735 * t703;
t860 = t646 * t650;
t647 = (t733 * t816 + t735 * t749) * t873 + t707 * t836 + t735 * t704;
t859 = t647 * t651;
t648 = (t733 * t815 + t735 * t751) * t872 + t708 * t836 + t735 * t705;
t858 = t648 * t652;
t801 = t739 * mrSges(3,1) - mrSges(3,2) * t737;
t661 = t801 * t736 - (t737 * mrSges(3,1) + t739 * mrSges(3,2)) * t834;
t857 = t661 * t774;
t822 = t736 * t738;
t682 = t733 * t740 + t735 * t822;
t833 = t734 * t739;
t663 = t737 * t682 + t735 * t833;
t856 = t663 * t717;
t855 = t663 * t721;
t820 = t736 * t747;
t687 = t733 * t753 + t735 * t820;
t826 = t734 * t752;
t667 = t746 * t687 + t735 * t826;
t854 = t667 * t718;
t853 = t667 * t722;
t819 = t736 * t749;
t688 = t733 * t755 + t735 * t819;
t825 = t734 * t754;
t668 = t748 * t688 + t735 * t825;
t852 = t668 * t719;
t851 = t668 * t723;
t818 = t736 * t751;
t689 = t733 * t757 + t735 * t818;
t824 = t734 * t756;
t669 = t750 * t689 + t735 * t824;
t850 = t669 * t720;
t849 = t669 * t724;
t800 = t752 * mrSges(3,1) - mrSges(3,2) * t746;
t670 = t800 * t736 - (t746 * mrSges(3,1) + t752 * mrSges(3,2)) * t831;
t848 = t670 * t774;
t799 = t754 * mrSges(3,1) - mrSges(3,2) * t748;
t671 = t799 * t736 - (t748 * mrSges(3,1) + t754 * mrSges(3,2)) * t829;
t847 = t671 * t774;
t798 = t756 * mrSges(3,1) - mrSges(3,2) * t750;
t672 = t798 * t736 - (t750 * mrSges(3,1) + t756 * mrSges(3,2)) * t827;
t846 = t672 * t774;
t710 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t845 = ((t801 + t870) * t740 + t710 * t738) * t734;
t844 = ((t800 + t870) * t753 + t710 * t747) * t734;
t843 = ((t799 + t870) * t755 + t710 * t749) * t734;
t842 = ((t798 + t870) * t757 + t710 * t751) * t734;
t712 = -mrSges(3,2) * pkin(6) + Ifges(3,6);
t713 = mrSges(3,1) * pkin(6) - Ifges(3,5);
t683 = t712 * t739 - t737 * t713;
t841 = t683 * t774;
t690 = t712 * t752 - t746 * t713;
t840 = t690 * t774;
t691 = t712 * t754 - t748 * t713;
t839 = t691 * t774;
t692 = t712 * t756 - t750 * t713;
t838 = t692 * t774;
t837 = t733 * t734;
t835 = t734 * t737;
t832 = t734 * t746;
t830 = t734 * t748;
t828 = t734 * t750;
t810 = -0.2e1 * pkin(2) * mrSges(3,2);
t809 = t717 * t869;
t808 = t721 * t869;
t807 = t718 * t866;
t806 = t722 * t866;
t805 = t719 * t864;
t804 = t723 * t864;
t803 = t720 * t862;
t802 = t724 * t862;
t764 = mrSges(4,2);
t765 = mrSges(4,1);
t797 = -t725 * t764 + t726 * t765;
t796 = 0.2e1 * pkin(6) * mrSges(3,3) + Ifges(3,1) + Ifges(2,3) + (pkin(2) ^ 2 + pkin(6) ^ 2) * m(3);
t791 = -t725 * t765 - t726 * t764;
t790 = pkin(3) * t835 - t701 * t736;
t789 = pkin(3) * t832 - t703 * t736;
t788 = pkin(3) * t830 - t704 * t736;
t787 = pkin(3) * t828 - t705 * t736;
t786 = Ifges(3,3) * t868 - t663 * t683;
t785 = Ifges(3,3) * t865 - t667 * t690;
t784 = Ifges(3,3) * t863 - t668 * t691;
t783 = Ifges(3,3) * t861 - t669 * t692;
t741 = -Ifges(3,1) + Ifges(3,2);
t758 = mrSges(3,1) * pkin(2);
t657 = t741 * t727 + 0.2e1 * (Ifges(3,4) * t737 + t758) * t739 + t737 * t810 + t796;
t782 = t641 * t841 - t657 * t663;
t658 = t741 * t730 + 0.2e1 * (Ifges(3,4) * t746 + t758) * t752 + t746 * t810 + t796;
t781 = t643 * t840 - t658 * t667;
t659 = t741 * t731 + 0.2e1 * (Ifges(3,4) * t748 + t758) * t754 + t748 * t810 + t796;
t780 = t644 * t839 - t659 * t668;
t660 = t741 * t732 + 0.2e1 * (Ifges(3,4) * t750 + t758) * t756 + t750 * t810 + t796;
t779 = t645 * t838 - t660 * t669;
t778 = t641 * t857 - t663 * t845;
t777 = t643 * t848 - t667 * t844;
t776 = t644 * t847 - t668 * t843;
t775 = t645 * t846 - t669 * t842;
t728 = m(1) + m(2) + m(3);
t686 = t733 * t818 - t735 * t757;
t685 = t733 * t819 - t735 * t755;
t684 = t733 * t820 - t735 * t753;
t681 = t733 * t822 - t735 * t740;
t666 = t750 * t686 + t733 * t824;
t665 = t748 * t685 + t733 * t825;
t664 = t746 * t684 + t733 * t826;
t662 = t737 * t681 + t733 * t833;
t656 = t708 * t735 + t787 * t733;
t655 = t707 * t735 + t788 * t733;
t654 = t706 * t735 + t789 * t733;
t653 = t702 * t735 + t790 * t733;
t640 = -t689 * t876 - t708 * t733 * t756 + (pkin(2) * t828 + t787 * t756) * t735;
t639 = -t688 * t877 - t707 * t733 * t754 + (pkin(2) * t830 + t788 * t754) * t735;
t638 = -t687 * t878 - t706 * t733 * t752 + (pkin(2) * t832 + t789 * t752) * t735;
t637 = -t682 * t879 - t702 * t733 * t739 + (pkin(2) * t835 + t790 * t739) * t735;
t636 = -(t686 * t724 - t720 * t827) * t876 + (t656 * t724 + t680 * t720) * t756 + (t736 * t720 + t724 * t837) * t880;
t635 = (t686 * t720 + t724 * t827) * t876 + (-t656 * t720 + t680 * t724) * t756 + (-t720 * t837 + t736 * t724) * t880;
t634 = -(t685 * t723 - t719 * t829) * t877 + (t655 * t723 + t679 * t719) * t754 + (t736 * t719 + t723 * t837) * t881;
t633 = (t685 * t719 + t723 * t829) * t877 + (-t655 * t719 + t679 * t723) * t754 + (-t719 * t837 + t736 * t723) * t881;
t632 = -(t684 * t722 - t718 * t831) * t878 + (t654 * t722 + t678 * t718) * t752 + (t736 * t718 + t722 * t837) * t882;
t631 = (t684 * t718 + t722 * t831) * t878 + (-t654 * t718 + t678 * t722) * t752 + (-t718 * t837 + t736 * t722) * t882;
t630 = -(t681 * t721 - t717 * t834) * t879 + (t653 * t721 + t677 * t717) * t739 + (t736 * t717 + t721 * t837) * t883;
t629 = (t681 * t717 + t721 * t834) * t879 + (-t653 * t717 + t677 * t721) * t739 + (-t717 * t837 + t736 * t721) * t883;
t628 = t669 * t884;
t627 = t668 * t885;
t626 = t667 * t886;
t625 = t663 * t887;
t624 = t861 * t884;
t623 = t863 * t885;
t622 = t865 * t886;
t621 = t868 * t887;
t620 = (t640 * t672 + t648 * t871 + t666 * t692) * t652;
t619 = (t639 * t671 + t647 * t871 + t665 * t691) * t651;
t618 = (t638 * t670 + t646 * t871 + t664 * t690) * t650;
t617 = (t637 * t661 + t642 * t871 + t662 * t683) * t649;
t616 = (t640 * t728 + t648 * t846 + t666 * t842) * t652;
t615 = (t639 * t728 + t647 * t847 + t665 * t843) * t651;
t614 = (t638 * t728 + t646 * t848 + t664 * t844) * t650;
t613 = (t637 * t728 + t642 * t857 + t662 * t845) * t649;
t612 = (t640 * t842 + t648 * t838 + t660 * t666) * t652;
t611 = (t639 * t843 + t647 * t839 + t659 * t665) * t651;
t610 = (t638 * t844 + t646 * t840 + t658 * t664) * t650;
t609 = (t637 * t845 + t642 * t841 + t657 * t662) * t649;
t608 = (t635 * t700 + t636 * t696) * t652;
t607 = (t633 * t699 + t634 * t695) * t651;
t606 = (t631 * t698 + t632 * t694) * t650;
t605 = (t629 * t697 + t630 * t693) * t649;
t604 = (t636 * t672 + t783 * t724) * t652;
t603 = (t635 * t672 - t783 * t720) * t652;
t602 = (t634 * t671 + t784 * t723) * t651;
t601 = (t633 * t671 - t784 * t719) * t651;
t600 = (t632 * t670 + t785 * t722) * t650;
t599 = (t631 * t670 - t785 * t718) * t650;
t598 = (t630 * t661 + t786 * t721) * t649;
t597 = (t629 * t661 - t786 * t717) * t649;
t596 = (t636 * t728 + t775 * t724) * t652;
t595 = (t635 * t728 - t775 * t720) * t652;
t594 = (t634 * t728 + t776 * t723) * t651;
t593 = (t633 * t728 - t776 * t719) * t651;
t592 = (t632 * t728 + t777 * t722) * t650;
t591 = (t631 * t728 - t777 * t718) * t650;
t590 = (t630 * t728 + t778 * t721) * t649;
t589 = (t629 * t728 - t778 * t717) * t649;
t588 = (t636 * t842 + t779 * t724) * t652;
t587 = (t635 * t842 - t779 * t720) * t652;
t586 = (t634 * t843 + t780 * t723) * t651;
t585 = (t633 * t843 - t780 * t719) * t651;
t584 = (t632 * t844 + t781 * t722) * t650;
t583 = (t631 * t844 - t781 * t718) * t650;
t582 = (t630 * t845 + t782 * t721) * t649;
t581 = (t629 * t845 - t782 * t717) * t649;
t580 = t624 * Ifges(3,3) + t608 * t672 - t628 * t692;
t579 = t623 * Ifges(3,3) + t607 * t671 - t627 * t691;
t578 = t622 * Ifges(3,3) + t606 * t670 - t626 * t690;
t577 = t608 * t728 + t624 * t672 - t628 * t842;
t576 = t607 * t728 + t623 * t671 - t627 * t843;
t575 = t606 * t728 + t622 * t670 - t626 * t844;
t574 = t621 * Ifges(3,3) + t605 * t661 - t625 * t683;
t573 = t605 * t728 + t621 * t661 - t625 * t845;
t572 = t608 * t842 + t624 * t692 - t628 * t660;
t571 = t607 * t843 + t623 * t691 - t627 * t659;
t570 = t606 * t844 + t622 * t690 - t626 * t658;
t569 = t605 * t845 + t621 * t683 - t625 * t657;
t1 = [m(4) + (-t588 * t849 + t596 * t636) * t652 + (-t586 * t851 + t594 * t634) * t651 + (-t584 * t853 + t592 * t632) * t650 + (-t582 * t855 + t590 * t630) * t649 + (t598 * t808 + t600 * t806 + t602 * t804 + t604 * t802) * t774, (t588 * t850 + t596 * t635) * t652 + (t586 * t852 + t594 * t633) * t651 + (t584 * t854 + t592 * t631) * t650 + (t582 * t856 + t590 * t629) * t649 + (-t598 * t809 - t600 * t807 - t602 * t805 - t604 * t803) * t774, (t588 * t666 + t596 * t640) * t652 + (t586 * t665 + t594 * t639) * t651 + (t584 * t664 + t592 * t638) * t650 + (t582 * t662 + t590 * t637) * t649 + (t598 * t867 + t600 * t860 + t602 * t859 + t604 * t858) * t774, -t582 * t625 - t584 * t626 - t586 * t627 - t588 * t628 + t590 * t605 + t592 * t606 + t594 * t607 + t596 * t608 + t598 * t621 + t600 * t622 + t602 * t623 + t604 * t624 + t791; (-t587 * t849 + t595 * t636) * t652 + (-t585 * t851 + t593 * t634) * t651 + (-t583 * t853 + t591 * t632) * t650 + (-t581 * t855 + t589 * t630) * t649 + (t597 * t808 + t599 * t806 + t601 * t804 + t603 * t802) * t774, m(4) + (t587 * t850 + t595 * t635) * t652 + (t585 * t852 + t593 * t633) * t651 + (t583 * t854 + t591 * t631) * t650 + (t581 * t856 + t589 * t629) * t649 + (-t597 * t809 - t599 * t807 - t601 * t805 - t603 * t803) * t774, (t587 * t666 + t595 * t640) * t652 + (t585 * t665 + t593 * t639) * t651 + (t583 * t664 + t591 * t638) * t650 + (t581 * t662 + t589 * t637) * t649 + (t597 * t867 + t599 * t860 + t601 * t859 + t603 * t858) * t774, -t581 * t625 - t583 * t626 - t585 * t627 - t587 * t628 + t589 * t605 + t591 * t606 + t593 * t607 + t595 * t608 + t597 * t621 + t599 * t622 + t601 * t623 + t603 * t624 + t797; (-t612 * t849 + t616 * t636) * t652 + (-t611 * t851 + t615 * t634) * t651 + (-t610 * t853 + t614 * t632) * t650 + (-t609 * t855 + t613 * t630) * t649 + (t617 * t808 + t618 * t806 + t619 * t804 + t620 * t802) * t774, (t612 * t850 + t616 * t635) * t652 + (t611 * t852 + t615 * t633) * t651 + (t610 * t854 + t614 * t631) * t650 + (t609 * t856 + t613 * t629) * t649 + (-t617 * t809 - t618 * t807 - t619 * t805 - t620 * t803) * t774, m(4) + (t612 * t666 + t616 * t640) * t652 + (t611 * t665 + t615 * t639) * t651 + (t610 * t664 + t614 * t638) * t650 + (t609 * t662 + t613 * t637) * t649 + (t617 * t867 + t618 * t860 + t619 * t859 + t620 * t858) * t774, t613 * t605 + t614 * t606 + t615 * t607 + t616 * t608 - t609 * t625 - t610 * t626 - t611 * t627 - t612 * t628 + t617 * t621 + t618 * t622 + t619 * t623 + t620 * t624; (-t572 * t849 + t577 * t636) * t652 + (-t571 * t851 + t576 * t634) * t651 + (-t570 * t853 + t575 * t632) * t650 + (-t569 * t855 + t573 * t630) * t649 + (t574 * t808 + t578 * t806 + t579 * t804 + t580 * t802) * t774 + t791, (t572 * t850 + t577 * t635) * t652 + (t571 * t852 + t576 * t633) * t651 + (t570 * t854 + t575 * t631) * t650 + (t569 * t856 + t573 * t629) * t649 + (-t574 * t809 - t578 * t807 - t579 * t805 - t580 * t803) * t774 + t797, (t572 * t666 + t577 * t640) * t652 + (t571 * t665 + t576 * t639) * t651 + (t570 * t664 + t575 * t638) * t650 + (t569 * t662 + t573 * t637) * t649 + (t574 * t867 + t578 * t860 + t579 * t859 + t580 * t858) * t774, -t569 * t625 - t570 * t626 - t571 * t627 - t572 * t628 + t573 * t605 + t574 * t621 + t575 * t606 + t576 * t607 + t577 * t608 + t578 * t622 + t579 * t623 + t580 * t624 + Ifges(4,3);];
MX  = t1;
