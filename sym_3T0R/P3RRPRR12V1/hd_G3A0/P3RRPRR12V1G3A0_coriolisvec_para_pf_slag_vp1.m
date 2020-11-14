% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRPRR12V1G3A0
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
% Datum: 2020-08-06 19:11
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR12V1G3A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G3A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V1G3A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G3A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G3A0_coriolisvec_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V1G3A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR12V1G3A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR12V1G3A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G3A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G3A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:09:49
% EndTime: 2020-08-06 19:09:53
% DurationCPUTime: 4.18s
% Computational Cost: add. (27468->300), mult. (41535->488), div. (6471->6), fcn. (34602->18), ass. (0->228)
t727 = cos(qJ(1,3));
t720 = sin(qJ(2,3));
t824 = qJ(3,3) * t720;
t721 = sin(qJ(1,3));
t837 = pkin(4) * t721;
t679 = t727 * t824 - t837;
t717 = legFrame(3,2);
t695 = sin(t717);
t698 = cos(t717);
t726 = cos(qJ(2,3));
t710 = t726 ^ 2;
t736 = pkin(1) + pkin(2);
t784 = t727 * t736;
t791 = t720 * t736;
t798 = t695 * qJ(3,3);
t639 = (t698 * t784 - t798) * t710 + (t679 * t698 + t695 * t791) * t726 + t798;
t795 = t698 * qJ(3,3);
t642 = (-t695 * t784 - t795) * t710 + (-t679 * t695 + t698 * t791) * t726 + t795;
t786 = t726 * t736;
t688 = t786 + t824;
t666 = pkin(4) * t727 + t688 * t721;
t733 = xDP(3);
t734 = xDP(2);
t735 = xDP(1);
t682 = 0.1e1 / t688;
t738 = 0.1e1 / qJ(3,3);
t801 = t682 * t738;
t624 = (-t666 * t726 * t733 + t639 * t735 + t642 * t734) * t801;
t861 = t624 * t736;
t729 = cos(qJ(1,2));
t722 = sin(qJ(2,2));
t826 = qJ(3,2) * t722;
t723 = sin(qJ(1,2));
t836 = pkin(4) * t723;
t680 = t729 * t826 - t836;
t718 = legFrame(2,2);
t696 = sin(t718);
t699 = cos(t718);
t728 = cos(qJ(2,2));
t711 = t728 ^ 2;
t781 = t729 * t736;
t789 = t722 * t736;
t797 = t696 * qJ(3,2);
t640 = (t699 * t781 - t797) * t711 + (t680 * t699 + t696 * t789) * t728 + t797;
t794 = t699 * qJ(3,2);
t643 = (-t696 * t781 - t794) * t711 + (-t680 * t696 + t699 * t789) * t728 + t794;
t783 = t728 * t736;
t689 = t783 + t826;
t667 = pkin(4) * t729 + t689 * t723;
t683 = 0.1e1 / t689;
t740 = 0.1e1 / qJ(3,2);
t800 = t683 * t740;
t625 = (-t667 * t728 * t733 + t640 * t735 + t643 * t734) * t800;
t860 = t625 * t736;
t731 = cos(qJ(1,1));
t724 = sin(qJ(2,1));
t828 = qJ(3,1) * t724;
t725 = sin(qJ(1,1));
t835 = pkin(4) * t725;
t681 = t731 * t828 - t835;
t719 = legFrame(1,2);
t697 = sin(t719);
t700 = cos(t719);
t730 = cos(qJ(2,1));
t712 = t730 ^ 2;
t778 = t731 * t736;
t787 = t724 * t736;
t796 = t697 * qJ(3,1);
t641 = (t700 * t778 - t796) * t712 + (t681 * t700 + t697 * t787) * t730 + t796;
t793 = t700 * qJ(3,1);
t644 = (-t697 * t778 - t793) * t712 + (-t681 * t697 + t700 * t787) * t730 + t793;
t780 = t730 * t736;
t690 = t780 + t828;
t668 = pkin(4) * t731 + t690 * t725;
t684 = 0.1e1 / t690;
t742 = 0.1e1 / qJ(3,1);
t799 = t684 * t742;
t626 = (-t668 * t730 * t733 + t641 * t735 + t644 * t734) * t799;
t859 = t626 * t736;
t823 = qJ(3,3) * t726;
t752 = -t791 + t823;
t858 = t752 * t624;
t825 = qJ(3,2) * t728;
t753 = -t789 + t825;
t857 = t753 * t625;
t827 = qJ(3,1) * t730;
t754 = -t787 + t827;
t856 = t754 * t626;
t855 = t688 * t727 - t837;
t854 = t689 * t729 - t836;
t853 = t690 * t731 - t835;
t636 = (-t727 * t733 + (t695 * t734 - t698 * t735) * t721) * t682;
t737 = qJ(3,3) ^ 2;
t748 = pkin(4) ^ 2;
t804 = t666 * t738;
t807 = (-t855 * t695 - t752 * t698) * t738;
t810 = (-t752 * t695 + t855 * t698) * t738;
t627 = -t733 * t804 + t734 * t807 + t735 * t810;
t603 = -t627 + t861;
t822 = t603 * t720;
t597 = (-t624 * t823 + t822) * pkin(4) + ((qJ(3,3) + t736) * (-qJ(3,3) + t736) * t710 + 0.2e1 * t786 * t824 + t737 + t748) * t636;
t817 = t636 * t720;
t630 = pkin(4) * t817;
t609 = t630 + t861;
t785 = t726 * t738;
t818 = t636 * t710;
t819 = t636 * t682;
t585 = -t597 * t785 * t819 + ((-(t630 + t603) * t786 + (pkin(4) * t818 - t822) * qJ(3,3)) * t624 + (t609 * t726 + t624 * t824) * t627) * t801;
t750 = (pkin(1) ^ 2);
t761 = -t750 + (-2 * pkin(1) - pkin(2)) * pkin(2);
t588 = (((-t737 + t761) * t624 + t736 * t627) * t624 + t609 * t627 + (pkin(4) * t858 - t597) * t636) * t738;
t594 = (-pkin(4) * t636 + 0.2e1 * t720 * t627 + 0.2e1 * t858) * t819;
t714 = rSges(3,3) + qJ(3,3);
t763 = -m(2) * rSges(2,1) * rSges(2,2) + Icges(2,4) - Icges(3,5);
t732 = pkin(1) + rSges(3,1);
t838 = m(3) * t732;
t672 = t714 * t838 + t763;
t612 = t672 * t624;
t621 = t624 ^ 2;
t842 = m(2) * rSges(2,3);
t762 = -rSges(2,2) * t842 + Icges(2,6) - Icges(3,6);
t841 = m(3) * t714;
t675 = rSges(3,2) * t841 + t762;
t678 = rSges(2,1) * t842 + rSges(3,2) * t838 - Icges(3,4) - Icges(2,5);
t657 = -t675 * t726 + t678 * t720;
t745 = rSges(2,2) ^ 2;
t747 = rSges(2,1) ^ 2;
t777 = -Icges(2,1) - Icges(3,1);
t755 = Icges(2,2) + Icges(3,3) + (-t745 + t747) * m(2) + t777;
t767 = rSges(3,3) + t732;
t768 = rSges(3,3) - t732;
t663 = -(qJ(3,3) + t767) * (qJ(3,3) + t768) * m(3) + t755;
t751 = -(rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - (rSges(2,3) ^ 2 + t745) * m(2) - Icges(1,3) + t777;
t766 = t627 * t841;
t769 = -0.2e1 * rSges(3,2) * m(3);
t743 = rSges(3,3) ^ 2;
t770 = rSges(3,2) ^ 2 + t743;
t792 = t720 * t726;
t831 = rSges(3,2) * t720;
t845 = -0.2e1 * t672;
t846 = 0.2e1 * rSges(3,3);
t849 = qJ(3,3) * t846 + t737;
t776 = (-t663 * t710 + t792 * t845 - (t770 + t849) * m(3) + t751) * t594 + t657 * t585 - m(3) * t588 * t831 - 0.4e1 * (-t612 + t766 / 0.2e1) * t818 + (-0.2e1 * (t624 * t663 - t627 * t838) * t817 - t624 * (t624 * t678 + t627 * t769)) * t726 - t621 * t675 * t720 + 0.2e1 * (-t612 + t766) * t636;
t852 = t721 * t776;
t637 = (-t729 * t733 + (t696 * t734 - t699 * t735) * t723) * t683;
t739 = qJ(3,2) ^ 2;
t803 = t667 * t740;
t806 = (-t854 * t696 - t753 * t699) * t740;
t809 = (-t753 * t696 + t854 * t699) * t740;
t628 = -t733 * t803 + t734 * t806 + t735 * t809;
t604 = -t628 + t860;
t821 = t604 * t722;
t598 = (-t625 * t825 + t821) * pkin(4) + ((qJ(3,2) + t736) * (-qJ(3,2) + t736) * t711 + 0.2e1 * t783 * t826 + t739 + t748) * t637;
t814 = t637 * t722;
t631 = pkin(4) * t814;
t610 = t631 + t860;
t782 = t728 * t740;
t815 = t637 * t711;
t816 = t637 * t683;
t586 = -t598 * t782 * t816 + ((-(t631 + t604) * t783 + (pkin(4) * t815 - t821) * qJ(3,2)) * t625 + (t610 * t728 + t625 * t826) * t628) * t800;
t589 = (((-t739 + t761) * t625 + t736 * t628) * t625 + t610 * t628 + (pkin(4) * t857 - t598) * t637) * t740;
t595 = (-pkin(4) * t637 + 0.2e1 * t722 * t628 + 0.2e1 * t857) * t816;
t715 = rSges(3,3) + qJ(3,2);
t673 = t715 * t838 + t763;
t613 = t673 * t625;
t622 = t625 ^ 2;
t840 = m(3) * t715;
t676 = rSges(3,2) * t840 + t762;
t658 = -t676 * t728 + t678 * t722;
t664 = -(qJ(3,2) + t767) * (qJ(3,2) + t768) * m(3) + t755;
t765 = t628 * t840;
t790 = t722 * t728;
t830 = rSges(3,2) * t722;
t844 = -0.2e1 * t673;
t848 = qJ(3,2) * t846 + t739;
t775 = (-t664 * t711 + t790 * t844 - (t770 + t848) * m(3) + t751) * t595 + t658 * t586 - m(3) * t589 * t830 - 0.4e1 * (-t613 + t765 / 0.2e1) * t815 + (-0.2e1 * (t625 * t664 - t628 * t838) * t814 - t625 * (t625 * t678 + t628 * t769)) * t728 - t622 * t676 * t722 + 0.2e1 * (-t613 + t765) * t637;
t851 = t723 * t775;
t638 = (-t731 * t733 + (t697 * t734 - t700 * t735) * t725) * t684;
t741 = qJ(3,1) ^ 2;
t802 = t668 * t742;
t805 = (-t853 * t697 - t754 * t700) * t742;
t808 = (-t754 * t697 + t853 * t700) * t742;
t629 = -t733 * t802 + t734 * t805 + t735 * t808;
t605 = -t629 + t859;
t820 = t605 * t724;
t599 = (-t626 * t827 + t820) * pkin(4) + ((qJ(3,1) + t736) * (-qJ(3,1) + t736) * t712 + 0.2e1 * t780 * t828 + t741 + t748) * t638;
t811 = t638 * t724;
t632 = pkin(4) * t811;
t611 = t632 + t859;
t779 = t730 * t742;
t812 = t638 * t712;
t813 = t638 * t684;
t587 = -t599 * t779 * t813 + ((-(t632 + t605) * t780 + (pkin(4) * t812 - t820) * qJ(3,1)) * t626 + (t611 * t730 + t626 * t828) * t629) * t799;
t590 = (((-t741 + t761) * t626 + t736 * t629) * t626 + t611 * t629 + (pkin(4) * t856 - t599) * t638) * t742;
t596 = (-pkin(4) * t638 + 0.2e1 * t724 * t629 + 0.2e1 * t856) * t813;
t716 = rSges(3,3) + qJ(3,1);
t674 = t716 * t838 + t763;
t614 = t674 * t626;
t623 = t626 ^ 2;
t839 = m(3) * t716;
t677 = rSges(3,2) * t839 + t762;
t659 = -t677 * t730 + t678 * t724;
t665 = -(qJ(3,1) + t767) * (qJ(3,1) + t768) * m(3) + t755;
t764 = t629 * t839;
t788 = t724 * t730;
t829 = rSges(3,2) * t724;
t843 = -0.2e1 * t674;
t847 = qJ(3,1) * t846 + t741;
t774 = (-t665 * t712 + t788 * t843 - (t770 + t847) * m(3) + t751) * t596 + t659 * t587 - m(3) * t590 * t829 - 0.4e1 * (-t614 + t764 / 0.2e1) * t812 + (-0.2e1 * (t626 * t665 - t629 * t838) * t811 - t626 * (t626 * t678 + t629 * t769)) * t730 - t623 * t677 * t724 + 0.2e1 * (-t614 + t764) * t638;
t850 = t725 * t774;
t633 = t636 ^ 2;
t756 = -(t745 + t747) * m(2) - Icges(3,2) - Icges(2,3);
t757 = t743 + t750 + ((2 * pkin(1) + rSges(3,1)) * rSges(3,1));
t773 = t657 * t594 + (-(t757 + t849) * m(3) + t756) * t585 + t588 * t838 + 0.2e1 * t624 * t766 + (t663 * t792 + t710 * t845 + t672) * t633;
t634 = t637 ^ 2;
t772 = t658 * t595 + (-(t757 + t848) * m(3) + t756) * t586 + t589 * t838 + 0.2e1 * t625 * t765 + (t664 * t790 + t711 * t844 + t673) * t634;
t635 = t638 ^ 2;
t771 = t659 * t596 + (-(t757 + t847) * m(3) + t756) * t587 + t590 * t838 + 0.2e1 * t626 * t764 + (t665 * t788 + t712 * t843 + t674) * t635;
t760 = t773 * t738;
t759 = t772 * t740;
t758 = t771 * t742;
t608 = t623 * t716 + (-t712 * t716 + t732 * t788 + t716) * t635;
t607 = t622 * t715 + (-t711 * t715 + t732 * t790 + t715) * t634;
t606 = t621 * t714 + (-t710 * t714 + t732 * t792 + t714) * t633;
t584 = (t587 * t732 - t596 * t829 - t590) * m(3);
t583 = (t586 * t732 - t595 * t830 - t589) * m(3);
t582 = (t585 * t732 - t594 * t831 - t588) * m(3);
t1 = [t582 * t810 + t583 * t809 + t584 * t808 + (-t606 * t810 - t607 * t809 - t608 * t808) * m(3) + (t641 * t758 - t700 * t850) * t684 + (t640 * t759 - t699 * t851) * t683 + (t639 * t760 - t698 * t852) * t682; t582 * t807 + t583 * t806 + t584 * t805 + (-t606 * t807 - t607 * t806 - t608 * t805) * m(3) + (t644 * t758 + t697 * t850) * t684 + (t643 * t759 + t696 * t851) * t683 + (t642 * t760 + t695 * t852) * t682; -t582 * t804 - t583 * t803 - t584 * t802 + (t606 * t804 + t607 * t803 + t608 * t802) * m(3) + (-t771 * t668 * t779 - t774 * t731) * t684 + (-t772 * t667 * t782 - t775 * t729) * t683 + (-t773 * t666 * t785 - t776 * t727) * t682;];
taucX  = t1;
