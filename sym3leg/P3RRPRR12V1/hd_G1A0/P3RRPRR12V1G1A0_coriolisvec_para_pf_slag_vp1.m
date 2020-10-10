% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRPRR12V1G1A0
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
% Datum: 2020-08-06 19:02
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR12V1G1A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G1A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V1G1A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G1A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G1A0_coriolisvec_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V1G1A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR12V1G1A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR12V1G1A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G1A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G1A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:01:33
% EndTime: 2020-08-06 19:01:37
% DurationCPUTime: 4.46s
% Computational Cost: add. (20109->294), mult. (32307->478), div. (4653->6), fcn. (29763->18), ass. (0->228)
t741 = sin(qJ(2,3));
t754 = xDP(3);
t755 = xDP(2);
t756 = xDP(1);
t759 = 0.1e1 / qJ(3,3);
t747 = cos(qJ(2,3));
t757 = pkin(1) + pkin(2);
t819 = t757 * t741;
t776 = t747 * qJ(3,3) - t819;
t816 = t757 * t747;
t825 = t741 * qJ(3,3);
t700 = t816 + t825;
t748 = cos(qJ(1,3));
t716 = t748 * pkin(4);
t742 = sin(qJ(1,3));
t666 = t700 * t742 + t716;
t846 = t742 * pkin(4);
t669 = t700 * t748 - t846;
t738 = legFrame(3,3);
t704 = sin(t738);
t707 = cos(t738);
t657 = t666 * t707 + t704 * t669;
t831 = t657 * t759;
t654 = -t704 * t666 + t707 * t669;
t834 = t654 * t759;
t630 = -t754 * t759 * t776 + t755 * t831 + t756 * t834;
t672 = -t704 * t742 + t707 * t748;
t688 = t742 * t825 + t716;
t689 = t748 * t825 - t846;
t642 = t672 * t816 - t704 * t688 + t689 * t707;
t673 = t704 * t748 + t707 * t742;
t643 = t673 * t816 + t688 * t707 + t704 * t689;
t694 = 0.1e1 / t700;
t627 = (t741 * t754 + (t642 * t756 + t643 * t755) * t747 * t694) * t759;
t867 = t757 * t627;
t873 = t630 - t867;
t876 = t741 * t873;
t743 = sin(qJ(2,2));
t761 = 0.1e1 / qJ(3,2);
t749 = cos(qJ(2,2));
t818 = t757 * t743;
t777 = t749 * qJ(3,2) - t818;
t815 = t757 * t749;
t823 = t743 * qJ(3,2);
t701 = t815 + t823;
t750 = cos(qJ(1,2));
t717 = t750 * pkin(4);
t744 = sin(qJ(1,2));
t667 = t701 * t744 + t717;
t845 = t744 * pkin(4);
t670 = t701 * t750 - t845;
t739 = legFrame(2,3);
t705 = sin(t739);
t708 = cos(t739);
t658 = t667 * t708 + t705 * t670;
t830 = t658 * t761;
t655 = -t705 * t667 + t708 * t670;
t833 = t655 * t761;
t631 = -t754 * t761 * t777 + t755 * t830 + t756 * t833;
t674 = -t705 * t744 + t708 * t750;
t690 = t744 * t823 + t717;
t691 = t750 * t823 - t845;
t644 = t674 * t815 - t705 * t690 + t691 * t708;
t675 = t705 * t750 + t708 * t744;
t645 = t675 * t815 + t690 * t708 + t705 * t691;
t695 = 0.1e1 / t701;
t628 = (t743 * t754 + (t644 * t756 + t645 * t755) * t749 * t695) * t761;
t866 = t757 * t628;
t872 = t631 - t866;
t875 = t743 * t872;
t745 = sin(qJ(2,1));
t763 = 0.1e1 / qJ(3,1);
t751 = cos(qJ(2,1));
t817 = t757 * t745;
t778 = t751 * qJ(3,1) - t817;
t814 = t757 * t751;
t821 = t745 * qJ(3,1);
t702 = t814 + t821;
t752 = cos(qJ(1,1));
t718 = t752 * pkin(4);
t746 = sin(qJ(1,1));
t668 = t702 * t746 + t718;
t844 = t746 * pkin(4);
t671 = t702 * t752 - t844;
t740 = legFrame(1,3);
t706 = sin(t740);
t709 = cos(t740);
t659 = t668 * t709 + t706 * t671;
t829 = t659 * t763;
t656 = -t706 * t668 + t709 * t671;
t832 = t656 * t763;
t632 = -t754 * t763 * t778 + t755 * t829 + t756 * t832;
t676 = -t706 * t746 + t709 * t752;
t692 = t746 * t821 + t718;
t693 = t752 * t821 - t844;
t646 = t676 * t814 - t706 * t692 + t693 * t709;
t677 = t706 * t752 + t709 * t746;
t647 = t677 * t814 + t692 * t709 + t706 * t693;
t696 = 0.1e1 / t702;
t629 = (t745 * t754 + (t646 * t756 + t647 * t755) * t751 * t696) * t763;
t865 = t757 * t629;
t871 = t632 - t865;
t874 = t745 * t871;
t870 = 0.2e1 * t627;
t869 = 0.2e1 * t628;
t868 = 0.2e1 * t629;
t758 = qJ(3,3) ^ 2;
t858 = 2 * rSges(3,3);
t864 = qJ(3,3) * t858 + t758;
t760 = qJ(3,2) ^ 2;
t863 = qJ(3,2) * t858 + t760;
t762 = qJ(3,1) ^ 2;
t862 = qJ(3,1) * t858 + t762;
t639 = (t672 * t755 - t673 * t756) * t694;
t731 = t747 ^ 2;
t769 = pkin(4) ^ 2;
t861 = pkin(4) * t876 + (-t758 - t769 - (qJ(3,3) + t757) * (-qJ(3,3) + t757) * t731) * t639;
t640 = (t674 * t755 - t675 * t756) * t695;
t732 = t749 ^ 2;
t860 = pkin(4) * t875 + (-t760 - t769 - (qJ(3,2) + t757) * (-qJ(3,2) + t757) * t732) * t640;
t641 = (t676 * t755 - t677 * t756) * t696;
t733 = t751 ^ 2;
t859 = pkin(4) * t874 + (-t762 - t769 - (qJ(3,1) + t757) * (-qJ(3,1) + t757) * t733) * t641;
t735 = rSges(3,3) + qJ(3,3);
t790 = -m(2) * rSges(2,1) * rSges(2,2) + Icges(2,4) - Icges(3,5);
t753 = pkin(1) + rSges(3,1);
t850 = m(3) * t753;
t681 = t735 * t850 + t790;
t857 = -0.2e1 * t681;
t736 = rSges(3,3) + qJ(3,2);
t682 = t736 * t850 + t790;
t856 = -0.2e1 * t682;
t737 = rSges(3,3) + qJ(3,1);
t683 = t737 * t850 + t790;
t855 = -0.2e1 * t683;
t854 = m(2) * rSges(2,3);
t853 = m(3) * t735;
t852 = m(3) * t736;
t851 = m(3) * t737;
t849 = pkin(4) * t627;
t848 = pkin(4) * t628;
t847 = pkin(4) * t629;
t843 = rSges(3,2) * t741;
t842 = rSges(3,2) * t743;
t841 = rSges(3,2) * t745;
t840 = t639 * t731;
t839 = t639 * t741;
t838 = t640 * t732;
t837 = t640 * t743;
t836 = t641 * t733;
t835 = t641 * t745;
t824 = t741 * t747;
t822 = t743 * t749;
t820 = t745 * t751;
t633 = pkin(4) * t839;
t612 = t633 + t867;
t797 = -0.2e1 * qJ(3,3) * (t639 * t819 - t849 / 0.2e1);
t588 = ((t731 * t797 + t861 * t747) * t639 + (-(t633 - t873) * t816 + (pkin(4) * t840 + t876) * qJ(3,3)) * t627 + (t612 * t747 + t627 * t825) * t630) * t694 * t759;
t771 = (pkin(1) ^ 2);
t788 = -t771 + (-2 * pkin(1) - pkin(2)) * pkin(2);
t591 = (((-t758 + t788) * t627 + t757 * t630) * t627 + t612 * t630 + (t747 * t797 + t776 * t849 + t861) * t639) * t759;
t597 = (-pkin(4) * t639 + 0.2e1 * t741 * t630 + t776 * t870) * t639 * t694;
t585 = (t588 * t753 - t597 * t843 - t591) * m(3);
t813 = t759 * t585;
t634 = pkin(4) * t837;
t613 = t634 + t866;
t798 = -0.2e1 * qJ(3,2) * (t640 * t818 - t848 / 0.2e1);
t589 = ((t732 * t798 + t860 * t749) * t640 + (-(t634 - t872) * t815 + (pkin(4) * t838 + t875) * qJ(3,2)) * t628 + (t613 * t749 + t628 * t823) * t631) * t695 * t761;
t592 = (((-t760 + t788) * t628 + t757 * t631) * t628 + t613 * t631 + (t749 * t798 + t777 * t848 + t860) * t640) * t761;
t598 = (-pkin(4) * t640 + 0.2e1 * t743 * t631 + t777 * t869) * t640 * t695;
t586 = (t589 * t753 - t598 * t842 - t592) * m(3);
t812 = t761 * t586;
t635 = pkin(4) * t835;
t614 = t635 + t865;
t799 = -0.2e1 * qJ(3,1) * (t641 * t817 - t847 / 0.2e1);
t590 = ((t733 * t799 + t859 * t751) * t641 + (-(t635 - t871) * t814 + (pkin(4) * t836 + t874) * qJ(3,1)) * t629 + (t614 * t751 + t629 * t821) * t632) * t696 * t763;
t593 = (((-t762 + t788) * t629 + t757 * t632) * t629 + t614 * t632 + (t751 * t799 + t778 * t847 + t859) * t641) * t763;
t599 = (-pkin(4) * t641 + 0.2e1 * t745 * t632 + t778 * t868) * t641 * t696;
t587 = (t590 * t753 - t599 * t841 - t593) * m(3);
t811 = t763 * t587;
t810 = -Icges(2,1) - Icges(3,1);
t615 = t681 * t627;
t624 = t627 ^ 2;
t789 = -rSges(2,2) * t854 + Icges(2,6) - Icges(3,6);
t684 = rSges(3,2) * t853 + t789;
t687 = rSges(2,1) * t854 + rSges(3,2) * t850 - Icges(3,4) - Icges(2,5);
t660 = -t684 * t747 + t741 * t687;
t766 = rSges(2,2) ^ 2;
t768 = rSges(2,1) ^ 2;
t782 = Icges(2,2) + Icges(3,3) + (-t766 + t768) * m(2) + t810;
t800 = rSges(3,3) + t753;
t801 = rSges(3,3) - t753;
t663 = -(qJ(3,3) + t800) * (qJ(3,3) + t801) * m(3) + t782;
t775 = -(rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) - (rSges(2,3) ^ 2 + t766) * m(2) - Icges(1,3) + t810;
t796 = t630 * t853;
t802 = -0.2e1 * rSges(3,2) * m(3);
t764 = rSges(3,3) ^ 2;
t803 = rSges(3,2) ^ 2 + t764;
t809 = (-t663 * t731 + t824 * t857 - (t803 + t864) * m(3) + t775) * t597 + t660 * t588 - m(3) * t591 * t843 - 0.4e1 * (-t615 + t796 / 0.2e1) * t840 + (-0.2e1 * (t663 * t627 - t630 * t850) * t839 - t627 * (t687 * t627 + t630 * t802)) * t747 - t624 * t684 * t741 + 0.2e1 * (-t615 + t796) * t639;
t616 = t682 * t628;
t625 = t628 ^ 2;
t685 = rSges(3,2) * t852 + t789;
t661 = -t685 * t749 + t743 * t687;
t664 = -(qJ(3,2) + t800) * (qJ(3,2) + t801) * m(3) + t782;
t795 = t631 * t852;
t808 = (-t664 * t732 + t822 * t856 - (t803 + t863) * m(3) + t775) * t598 + t661 * t589 - m(3) * t592 * t842 - 0.4e1 * (-t616 + t795 / 0.2e1) * t838 + (-0.2e1 * (t664 * t628 - t631 * t850) * t837 - t628 * (t687 * t628 + t631 * t802)) * t749 - t625 * t685 * t743 + 0.2e1 * (-t616 + t795) * t640;
t617 = t683 * t629;
t626 = t629 ^ 2;
t686 = rSges(3,2) * t851 + t789;
t662 = -t686 * t751 + t745 * t687;
t665 = -(qJ(3,1) + t800) * (qJ(3,1) + t801) * m(3) + t782;
t794 = t632 * t851;
t807 = (-t665 * t733 + t820 * t855 - (t803 + t862) * m(3) + t775) * t599 + t662 * t590 - m(3) * t593 * t841 - 0.4e1 * (-t617 + t794 / 0.2e1) * t836 + (-0.2e1 * (t665 * t629 - t632 * t850) * t835 - t629 * (t687 * t629 + t632 * t802)) * t751 - t626 * t686 * t745 + 0.2e1 * (-t617 + t794) * t641;
t636 = t639 ^ 2;
t786 = -(t766 + t768) * m(2) - Icges(3,2) - Icges(2,3);
t787 = t764 + t771 + (2 * pkin(1) + rSges(3,1)) * rSges(3,1);
t806 = t660 * t597 + (-(t787 + t864) * m(3) + t786) * t588 + t591 * t850 + t796 * t870 + (t663 * t824 + t731 * t857 + t681) * t636;
t637 = t640 ^ 2;
t805 = t661 * t598 + (-(t787 + t863) * m(3) + t786) * t589 + t592 * t850 + t795 * t869 + (t664 * t822 + t732 * t856 + t682) * t637;
t638 = t641 ^ 2;
t804 = t662 * t599 + (-(t787 + t862) * m(3) + t786) * t590 + t593 * t850 + t794 * t868 + (t665 * t820 + t733 * t855 + t683) * t638;
t785 = t806 * t759 * t747;
t784 = t805 * t761 * t749;
t783 = t804 * t763 * t751;
t608 = t626 * t737 + (-t737 * t733 + t753 * t820 + t737) * t638;
t607 = t625 * t736 + (-t736 * t732 + t753 * t822 + t736) * t637;
t606 = t624 * t735 + (-t735 * t731 + t753 * t824 + t735) * t636;
t1 = [t654 * t813 + t655 * t812 + t656 * t811 + (-t606 * t834 - t607 * t833 - t608 * t832) * m(3) + (t646 * t783 - t807 * t677) * t696 + (t644 * t784 - t808 * t675) * t695 + (t642 * t785 - t809 * t673) * t694; t657 * t813 + t658 * t812 + t659 * t811 + (-t606 * t831 - t607 * t830 - t608 * t829) * m(3) + (t647 * t783 + t807 * t676) * t696 + (t645 * t784 + t808 * t674) * t695 + (t643 * t785 + t809 * t672) * t694; (t804 * t745 - (-m(3) * t608 + t587) * t778) * t763 + (t805 * t743 - (-m(3) * t607 + t586) * t777) * t761 + (t806 * t741 - (-m(3) * t606 + t585) * t776) * t759;];
taucX  = t1;
