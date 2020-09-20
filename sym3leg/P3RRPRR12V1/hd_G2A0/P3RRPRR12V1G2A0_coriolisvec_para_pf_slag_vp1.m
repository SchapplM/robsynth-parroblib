% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRPRR12V1G2A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4]';
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
% Datum: 2020-08-06 19:07
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR12V1G2A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G2A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V1G2A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G2A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G2A0_coriolisvec_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V1G2A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR12V1G2A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR12V1G2A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G2A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G2A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:05:32
% EndTime: 2020-08-06 19:05:36
% DurationCPUTime: 3.96s
% Computational Cost: add. (27468->300), mult. (42123->485), div. (6471->6), fcn. (35190->18), ass. (0->228)
t713 = sin(qJ(1,3));
t712 = sin(qJ(2,3));
t793 = t712 * qJ(3,3);
t719 = cos(qJ(1,3));
t826 = t719 * pkin(4);
t674 = t713 * t793 + t826;
t709 = legFrame(3,2);
t687 = sin(t709);
t690 = cos(t709);
t718 = cos(qJ(2,3));
t702 = t718 ^ 2;
t728 = pkin(1) + pkin(2);
t778 = t728 * t712;
t790 = t713 * t728;
t799 = t687 * qJ(3,3);
t631 = (t690 * t790 - t799) * t702 + (t674 * t690 + t687 * t778) * t718 + t799;
t796 = t690 * qJ(3,3);
t634 = (-t687 * t790 - t796) * t702 + (-t687 * t674 + t690 * t778) * t718 + t796;
t725 = xDP(3);
t726 = xDP(2);
t727 = xDP(1);
t775 = t728 * t718;
t683 = t775 + t793;
t677 = 0.1e1 / t683;
t730 = 0.1e1 / qJ(3,3);
t802 = t677 * t730;
t750 = -t713 * pkin(4) + t683 * t719;
t808 = t750 * t718;
t616 = (t631 * t727 + t634 * t726 + t725 * t808) * t802;
t847 = t616 * t728;
t715 = sin(qJ(1,2));
t714 = sin(qJ(2,2));
t789 = t714 * qJ(3,2);
t721 = cos(qJ(1,2));
t825 = t721 * pkin(4);
t675 = t715 * t789 + t825;
t710 = legFrame(2,2);
t688 = sin(t710);
t691 = cos(t710);
t720 = cos(qJ(2,2));
t703 = t720 ^ 2;
t777 = t728 * t714;
t786 = t715 * t728;
t798 = t688 * qJ(3,2);
t632 = (t691 * t786 - t798) * t703 + (t675 * t691 + t688 * t777) * t720 + t798;
t795 = t691 * qJ(3,2);
t635 = (-t688 * t786 - t795) * t703 + (-t688 * t675 + t691 * t777) * t720 + t795;
t774 = t728 * t720;
t684 = t774 + t789;
t678 = 0.1e1 / t684;
t732 = 0.1e1 / qJ(3,2);
t801 = t678 * t732;
t749 = -t715 * pkin(4) + t684 * t721;
t807 = t749 * t720;
t617 = (t632 * t727 + t635 * t726 + t725 * t807) * t801;
t846 = t617 * t728;
t717 = sin(qJ(1,1));
t716 = sin(qJ(2,1));
t785 = t716 * qJ(3,1);
t723 = cos(qJ(1,1));
t824 = t723 * pkin(4);
t676 = t717 * t785 + t824;
t711 = legFrame(1,2);
t689 = sin(t711);
t692 = cos(t711);
t722 = cos(qJ(2,1));
t704 = t722 ^ 2;
t776 = t728 * t716;
t782 = t717 * t728;
t797 = t689 * qJ(3,1);
t633 = (t692 * t782 - t797) * t704 + (t676 * t692 + t689 * t776) * t722 + t797;
t794 = t692 * qJ(3,1);
t636 = (-t689 * t782 - t794) * t704 + (-t689 * t676 + t692 * t776) * t722 + t794;
t773 = t728 * t722;
t685 = t773 + t785;
t679 = 0.1e1 / t685;
t734 = 0.1e1 / qJ(3,1);
t800 = t679 * t734;
t748 = -t717 * pkin(4) + t685 * t723;
t806 = t748 * t722;
t618 = (t633 * t727 + t636 * t726 + t725 * t806) * t800;
t845 = t618 * t728;
t781 = t718 * qJ(3,3);
t744 = -t778 + t781;
t844 = t744 * t616;
t780 = t720 * qJ(3,2);
t745 = -t777 + t780;
t843 = t745 * t617;
t779 = t722 * qJ(3,1);
t746 = -t776 + t779;
t842 = t746 * t618;
t628 = (-t713 * t725 + (-t687 * t726 + t690 * t727) * t719) * t677;
t729 = qJ(3,3) ^ 2;
t740 = pkin(4) ^ 2;
t805 = t750 * t730;
t658 = t683 * t713 + t826;
t644 = -t658 * t687 - t690 * t744;
t813 = t644 * t730;
t643 = t658 * t690 - t687 * t744;
t814 = t643 * t730;
t619 = t725 * t805 + t726 * t813 + t727 * t814;
t595 = -t619 + t847;
t792 = t712 * t595;
t589 = (-t616 * t781 + t792) * pkin(4) + ((qJ(3,3) + t728) * (-qJ(3,3) + t728) * t702 + 0.2e1 * t775 * t793 + t729 + t740) * t628;
t819 = t628 * t712;
t622 = pkin(4) * t819;
t601 = t622 + t847;
t820 = t628 * t702;
t577 = (-t589 * t718 * t628 + (-(t622 + t595) * t775 + (pkin(4) * t820 - t792) * qJ(3,3)) * t616 + (t601 * t718 + t616 * t793) * t619) * t802;
t742 = (pkin(1) ^ 2);
t756 = -t742 + (-2 * pkin(1) - pkin(2)) * pkin(2);
t580 = (((-t729 + t756) * t616 + t728 * t619) * t616 + t601 * t619 + (pkin(4) * t844 - t589) * t628) * t730;
t586 = (-pkin(4) * t628 + 0.2e1 * t712 * t619 + 0.2e1 * t844) * t628 * t677;
t706 = rSges(3,3) + qJ(3,3);
t758 = -m(2) * rSges(2,1) * rSges(2,2) + Icges(2,4) - Icges(3,5);
t724 = pkin(1) + rSges(3,1);
t827 = m(3) * t724;
t667 = t706 * t827 + t758;
t604 = t667 * t616;
t613 = t616 ^ 2;
t831 = m(2) * rSges(2,3);
t757 = -rSges(2,2) * t831 + Icges(2,6) - Icges(3,6);
t830 = m(3) * t706;
t670 = rSges(3,2) * t830 + t757;
t673 = rSges(2,1) * t831 + rSges(3,2) * t827 - Icges(3,4) - Icges(2,5);
t649 = -t670 * t718 + t673 * t712;
t737 = rSges(2,2) ^ 2;
t739 = rSges(2,1) ^ 2;
t769 = -Icges(2,1) - Icges(3,1);
t747 = Icges(2,2) + Icges(3,3) + (-t737 + t739) * m(2) + t769;
t762 = rSges(3,3) + t724;
t763 = rSges(3,3) - t724;
t655 = -(qJ(3,3) + t762) * (qJ(3,3) + t763) * m(3) + t747;
t743 = -(rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - (rSges(2,3) ^ 2 + t737) * m(2) - Icges(1,3) + t769;
t761 = t619 * t830;
t764 = -0.2e1 * rSges(3,2) * m(3);
t735 = rSges(3,3) ^ 2;
t765 = rSges(3,2) ^ 2 + t735;
t791 = t712 * t718;
t823 = rSges(3,2) * t712;
t834 = -0.2e1 * t667;
t835 = 0.2e1 * rSges(3,3);
t838 = qJ(3,3) * t835 + t729;
t768 = (-t655 * t702 + t791 * t834 - (t765 + t838) * m(3) + t743) * t586 + t649 * t577 - m(3) * t580 * t823 - 0.4e1 * (-t604 + t761 / 0.2e1) * t820 + (-0.2e1 * (t655 * t616 - t619 * t827) * t819 - (t673 * t616 + t619 * t764) * t616) * t718 - t670 * t613 * t712 + 0.2e1 * t628 * (-t604 + t761);
t841 = t719 * t768;
t629 = (-t715 * t725 + (-t688 * t726 + t691 * t727) * t721) * t678;
t731 = qJ(3,2) ^ 2;
t804 = t749 * t732;
t659 = t684 * t715 + t825;
t646 = -t659 * t688 - t691 * t745;
t811 = t646 * t732;
t645 = t659 * t691 - t688 * t745;
t812 = t645 * t732;
t620 = t725 * t804 + t726 * t811 + t727 * t812;
t596 = -t620 + t846;
t788 = t714 * t596;
t590 = (-t617 * t780 + t788) * pkin(4) + ((qJ(3,2) + t728) * (-qJ(3,2) + t728) * t703 + 0.2e1 * t774 * t789 + t731 + t740) * t629;
t817 = t629 * t714;
t623 = pkin(4) * t817;
t602 = t623 + t846;
t818 = t629 * t703;
t578 = (-t590 * t720 * t629 + (-(t623 + t596) * t774 + (pkin(4) * t818 - t788) * qJ(3,2)) * t617 + (t602 * t720 + t617 * t789) * t620) * t801;
t581 = (((-t731 + t756) * t617 + t728 * t620) * t617 + t602 * t620 + (pkin(4) * t843 - t590) * t629) * t732;
t587 = (-pkin(4) * t629 + 0.2e1 * t714 * t620 + 0.2e1 * t843) * t629 * t678;
t707 = rSges(3,3) + qJ(3,2);
t668 = t707 * t827 + t758;
t605 = t668 * t617;
t614 = t617 ^ 2;
t829 = m(3) * t707;
t671 = rSges(3,2) * t829 + t757;
t650 = -t671 * t720 + t673 * t714;
t656 = -(qJ(3,2) + t762) * (qJ(3,2) + t763) * m(3) + t747;
t760 = t620 * t829;
t787 = t714 * t720;
t822 = rSges(3,2) * t714;
t833 = -0.2e1 * t668;
t837 = qJ(3,2) * t835 + t731;
t767 = (-t656 * t703 + t787 * t833 - (t765 + t837) * m(3) + t743) * t587 + t650 * t578 - m(3) * t581 * t822 - 0.4e1 * (-t605 + t760 / 0.2e1) * t818 + (-0.2e1 * (t656 * t617 - t620 * t827) * t817 - (t673 * t617 + t620 * t764) * t617) * t720 - t671 * t614 * t714 + 0.2e1 * t629 * (-t605 + t760);
t840 = t721 * t767;
t630 = (-t717 * t725 + (-t689 * t726 + t692 * t727) * t723) * t679;
t733 = qJ(3,1) ^ 2;
t803 = t748 * t734;
t660 = t685 * t717 + t824;
t648 = -t660 * t689 - t692 * t746;
t809 = t648 * t734;
t647 = t660 * t692 - t689 * t746;
t810 = t647 * t734;
t621 = t725 * t803 + t726 * t809 + t727 * t810;
t597 = -t621 + t845;
t784 = t716 * t597;
t591 = (-t618 * t779 + t784) * pkin(4) + ((qJ(3,1) + t728) * (-qJ(3,1) + t728) * t704 + 0.2e1 * t773 * t785 + t733 + t740) * t630;
t815 = t630 * t716;
t624 = pkin(4) * t815;
t603 = t624 + t845;
t816 = t630 * t704;
t579 = (-t591 * t722 * t630 + (-(t624 + t597) * t773 + (pkin(4) * t816 - t784) * qJ(3,1)) * t618 + (t603 * t722 + t618 * t785) * t621) * t800;
t582 = (((-t733 + t756) * t618 + t728 * t621) * t618 + t603 * t621 + (pkin(4) * t842 - t591) * t630) * t734;
t588 = (-pkin(4) * t630 + 0.2e1 * t716 * t621 + 0.2e1 * t842) * t630 * t679;
t708 = rSges(3,3) + qJ(3,1);
t669 = t708 * t827 + t758;
t606 = t669 * t618;
t615 = t618 ^ 2;
t828 = m(3) * t708;
t672 = rSges(3,2) * t828 + t757;
t651 = -t672 * t722 + t673 * t716;
t657 = -(qJ(3,1) + t762) * (qJ(3,1) + t763) * m(3) + t747;
t759 = t621 * t828;
t783 = t716 * t722;
t821 = rSges(3,2) * t716;
t832 = -0.2e1 * t669;
t836 = qJ(3,1) * t835 + t733;
t766 = (-t657 * t704 + t783 * t832 - (t765 + t836) * m(3) + t743) * t588 + t651 * t579 - m(3) * t582 * t821 - 0.4e1 * (-t606 + t759 / 0.2e1) * t816 + (-0.2e1 * (t657 * t618 - t621 * t827) * t815 - (t673 * t618 + t621 * t764) * t618) * t722 - t672 * t615 * t716 + 0.2e1 * t630 * (-t606 + t759);
t839 = t723 * t766;
t772 = t730 * (t577 * t724 - t586 * t823 - t580) * m(3);
t771 = t732 * (t578 * t724 - t587 * t822 - t581) * m(3);
t770 = t734 * (t579 * t724 - t588 * t821 - t582) * m(3);
t625 = t628 ^ 2;
t751 = -(t737 + t739) * m(2) - Icges(3,2) - Icges(2,3);
t752 = t735 + t742 + ((2 * pkin(1) + rSges(3,1)) * rSges(3,1));
t755 = (t649 * t586 + (-(t752 + t838) * m(3) + t751) * t577 + t580 * t827 + 0.2e1 * t616 * t761 + (t655 * t791 + t702 * t834 + t667) * t625) * t730;
t626 = t629 ^ 2;
t754 = (t650 * t587 + (-(t752 + t837) * m(3) + t751) * t578 + t581 * t827 + 0.2e1 * t617 * t760 + (t656 * t787 + t703 * t833 + t668) * t626) * t732;
t627 = t630 ^ 2;
t753 = (t651 * t588 + (-(t752 + t836) * m(3) + t751) * t579 + t582 * t827 + 0.2e1 * t618 * t759 + (t657 * t783 + t704 * t832 + t669) * t627) * t734;
t600 = t615 * t708 + (-t708 * t704 + t724 * t783 + t708) * t627;
t599 = t614 * t707 + (-t707 * t703 + t724 * t787 + t707) * t626;
t598 = t613 * t706 + (-t706 * t702 + t724 * t791 + t706) * t625;
t1 = [t643 * t772 + t645 * t771 + t647 * t770 + (-t598 * t814 - t599 * t812 - t600 * t810) * m(3) + (t633 * t753 + t692 * t839) * t679 + (t632 * t754 + t691 * t840) * t678 + (t631 * t755 + t690 * t841) * t677; t644 * t772 + t646 * t771 + t648 * t770 + (-t598 * t813 - t599 * t811 - t600 * t809) * m(3) + (t636 * t753 - t689 * t839) * t679 + (t635 * t754 - t688 * t840) * t678 + (t634 * t755 - t687 * t841) * t677; t750 * t772 + t749 * t771 + t748 * t770 + (-t598 * t805 - t599 * t804 - t600 * t803) * m(3) + (-t766 * t717 + t753 * t806) * t679 + (-t767 * t715 + t754 * t807) * t678 + (-t768 * t713 + t755 * t808) * t677;];
taucX  = t1;
