% Calculate inertia matrix for parallel robot
% P3RPRRR6V1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:37
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPRRR6V1G2A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G2A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G2A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G2A0_inertia_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR6V1G2A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR6V1G2A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR6V1G2A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G2A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G2A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:35:49
% EndTime: 2020-08-06 18:35:51
% DurationCPUTime: 1.92s
% Computational Cost: add. (4206->263), mult. (6048->449), div. (531->10), fcn. (3690->56), ass. (0->186)
t681 = sin(pkin(7));
t702 = -pkin(6) - pkin(5);
t643 = t702 * t681 - pkin(1);
t687 = sin(qJ(1,3));
t627 = t643 * t687;
t682 = cos(pkin(7));
t693 = cos(qJ(1,3));
t750 = t693 * t702;
t751 = t693 * t681;
t796 = t627 - pkin(2) * t751 - (t687 * pkin(2) + t750) * t682;
t689 = sin(qJ(1,2));
t628 = t643 * t689;
t695 = cos(qJ(1,2));
t748 = t695 * t702;
t749 = t695 * t681;
t795 = t628 - pkin(2) * t749 - (t689 * pkin(2) + t748) * t682;
t691 = sin(qJ(1,1));
t629 = t643 * t691;
t697 = cos(qJ(1,1));
t746 = t697 * t702;
t747 = t697 * t681;
t794 = t629 - pkin(2) * t747 - (t691 * pkin(2) + t746) * t682;
t744 = 2 * m(3);
t793 = m(3) / 0.2e1;
t792 = Icges(3,2) / 0.2e1;
t690 = sin(qJ(3,1));
t696 = cos(qJ(3,1));
t638 = t696 * rSges(3,1) - t690 * rSges(3,2);
t688 = sin(qJ(3,2));
t694 = cos(qJ(3,2));
t637 = t694 * rSges(3,1) - t688 * rSges(3,2);
t686 = sin(qJ(3,3));
t692 = cos(qJ(3,3));
t636 = t692 * rSges(3,1) - t686 * rSges(3,2);
t791 = -0.2e1 * pkin(1);
t790 = -0.2e1 * pkin(2);
t789 = 0.2e1 * pkin(2);
t788 = 0.2e1 * t702;
t787 = m(3) * rSges(3,2);
t708 = rSges(3,2) ^ 2;
t709 = rSges(3,1) ^ 2;
t786 = (-t708 + t709) * t793 + t792 - Icges(3,1) / 0.2e1;
t785 = m(3) * t636;
t784 = m(3) * t637;
t783 = m(3) * t638;
t782 = pkin(1) * t681;
t659 = t682 * pkin(1);
t652 = t692 * pkin(3) + pkin(2);
t653 = t694 * pkin(3) + pkin(2);
t654 = t696 * pkin(3) + pkin(2);
t710 = 0.1e1 / pkin(3);
t760 = t652 * t687;
t775 = ((t750 + t760) * t682 - t627 + t652 * t751) * t710;
t759 = t653 * t689;
t774 = ((t748 + t759) * t682 - t628 + t653 * t749) * t710;
t758 = t654 * t691;
t773 = ((t746 + t758) * t682 - t629 + t654 * t747) * t710;
t669 = qJ(1,3) + pkin(7);
t656 = sin(t669);
t705 = 0.2e1 * qJ(3,3);
t672 = sin(t705);
t738 = pkin(7) + qJ(3,3);
t741 = -pkin(7) + qJ(3,3);
t772 = (t656 * t788 + cos(t669) * t790 + t693 * t791 + (-cos(qJ(1,3) - t741) - cos(qJ(1,3) + t738)) * pkin(3)) / (t686 * t789 + pkin(3) * t672 + (sin(t738) + sin(t741)) * pkin(1));
t670 = qJ(1,2) + pkin(7);
t657 = sin(t670);
t706 = 0.2e1 * qJ(3,2);
t673 = sin(t706);
t739 = pkin(7) + qJ(3,2);
t742 = -pkin(7) + qJ(3,2);
t771 = (t657 * t788 + cos(t670) * t790 + t695 * t791 + (-cos(qJ(1,2) - t742) - cos(qJ(1,2) + t739)) * pkin(3)) / (t688 * t789 + pkin(3) * t673 + (sin(t739) + sin(t742)) * pkin(1));
t671 = qJ(1,1) + pkin(7);
t658 = sin(t671);
t707 = 0.2e1 * qJ(3,1);
t674 = sin(t707);
t740 = pkin(7) + qJ(3,1);
t743 = -pkin(7) + qJ(3,1);
t770 = (t658 * t788 + cos(t671) * t790 + t697 * t791 + (-cos(qJ(1,1) - t743) - cos(qJ(1,1) + t740)) * pkin(3)) / (t690 * t789 + pkin(3) * t674 + (sin(t740) + sin(t743)) * pkin(1));
t683 = legFrame(3,2);
t645 = t683 + t669;
t646 = -t683 + t669;
t769 = -sin(t645) / 0.2e1 + sin(t646) / 0.2e1;
t684 = legFrame(2,2);
t647 = t684 + t670;
t648 = -t684 + t670;
t768 = -sin(t647) / 0.2e1 + sin(t648) / 0.2e1;
t685 = legFrame(1,2);
t649 = t685 + t671;
t650 = -t685 + t671;
t767 = -sin(t649) / 0.2e1 + sin(t650) / 0.2e1;
t766 = cos(t646) / 0.2e1 + cos(t645) / 0.2e1;
t765 = cos(t648) / 0.2e1 + cos(t647) / 0.2e1;
t764 = cos(t650) / 0.2e1 + cos(t649) / 0.2e1;
t630 = 0.1e1 / (t659 + t652);
t675 = 0.1e1 / t686;
t763 = t630 * t675;
t631 = 0.1e1 / (t659 + t653);
t676 = 0.1e1 / t688;
t762 = t631 * t676;
t632 = 0.1e1 / (t659 + t654);
t677 = 0.1e1 / t690;
t761 = t632 * t677;
t660 = sin(t683);
t757 = t660 * t686;
t661 = sin(t684);
t756 = t661 * t688;
t662 = sin(t685);
t755 = t662 * t690;
t663 = cos(t683);
t754 = t663 * t686;
t664 = cos(t684);
t753 = t664 * t688;
t665 = cos(t685);
t752 = t665 * t690;
t745 = t708 + t709;
t737 = (t687 * t682 + t751) * t692 ^ 2 * pkin(3);
t736 = (t689 * t682 + t749) * t694 ^ 2 * pkin(3);
t735 = (t691 * t682 + t747) * t696 ^ 2 * pkin(3);
t734 = t660 * t775;
t733 = t663 * t775;
t732 = ((t652 * t693 - t687 * t702) * t682 - t643 * t693 - t681 * t760) * t675 * t692;
t731 = t661 * t774;
t730 = t664 * t774;
t729 = ((t653 * t695 - t689 * t702) * t682 - t643 * t695 - t681 * t759) * t676 * t694;
t728 = t662 * t773;
t727 = t665 * t773;
t726 = ((t654 * t697 - t691 * t702) * t682 - t643 * t697 - t681 * t758) * t677 * t696;
t725 = t710 * t772;
t724 = t710 * t771;
t723 = t710 * t770;
t698 = rSges(3,3) + pkin(5);
t722 = t698 + t782;
t721 = t775 * t785;
t720 = t774 * t784;
t719 = t773 * t783;
t620 = -t722 * t787 + Icges(3,6);
t621 = t722 * rSges(3,1) * m(3) - Icges(3,5);
t608 = t620 * t692 - t686 * t621;
t718 = t608 * t675 * t775;
t609 = t620 * t694 - t688 * t621;
t717 = t609 * t676 * t774;
t610 = t620 * t696 - t690 * t621;
t716 = t610 * t677 * t773;
t711 = pkin(1) ^ 2;
t712 = Icges(1,3) + Icges(2,3) + (0.2e1 * pkin(2) ^ 2 + 0.2e1 * t698 ^ 2 + 0.2e1 * t711 + t745) * t793 + t698 * t782 * t744 + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + t792 + Icges(3,1) / 0.2e1 + (t711 + (-0.2e1 * t782 + rSges(2,2)) * rSges(2,2) + (0.2e1 * t659 + rSges(2,1)) * rSges(2,1)) * m(2);
t703 = m(2) + m(3);
t655 = -rSges(3,1) * t787 + Icges(3,4);
t651 = t659 + pkin(2);
t642 = t745 * m(3) + Icges(3,3);
t598 = t665 * t735 + (pkin(3) * t755 - t794 * t665) * t696 + t651 * t755;
t597 = -t662 * t735 + (pkin(3) * t752 + t794 * t662) * t696 + t651 * t752;
t596 = t664 * t736 + (pkin(3) * t756 - t795 * t664) * t694 + t651 * t756;
t595 = -t661 * t736 + (pkin(3) * t753 + t795 * t661) * t694 + t651 * t753;
t594 = t663 * t737 + (pkin(3) * t757 - t796 * t663) * t692 + t651 * t757;
t593 = -t660 * t737 + (pkin(3) * t754 + t796 * t660) * t692 + t651 * t754;
t592 = cos(t707) * t786 + t655 * t674 + (t638 * t659 + (t659 + t638) * pkin(2)) * t744 + t712;
t591 = cos(t706) * t786 + t655 * t673 + (t637 * t659 + (t659 + t637) * pkin(2)) * t744 + t712;
t590 = cos(t705) * t786 + t655 * t672 + (t636 * t659 + (t659 + t636) * pkin(2)) * t744 + t712;
t589 = t632 * t703 * t726 + t723 * t783;
t588 = t631 * t703 * t729 + t724 * t784;
t587 = t630 * t703 * t732 + t725 * t785;
t586 = (t598 * t703 - t665 * t719) * t761;
t585 = (t597 * t703 + t662 * t719) * t761;
t584 = (t596 * t703 - t664 * t720) * t762;
t583 = (t595 * t703 + t661 * t720) * t762;
t582 = (t594 * t703 - t663 * t721) * t763;
t581 = (t593 * t703 + t660 * t721) * t763;
t580 = t642 * t723 + (-t610 * t658 + t726 * t783) * t632;
t579 = t642 * t724 + (-t609 * t657 + t729 * t784) * t631;
t578 = t642 * t725 + (-t608 * t656 + t732 * t785) * t630;
t577 = -t658 * t632 * t592 + t610 * t723;
t576 = -t657 * t631 * t591 + t609 * t724;
t575 = -t656 * t630 * t590 + t608 * t725;
t574 = (t592 * t764 - t665 * t716) * t632;
t573 = (t591 * t765 - t664 * t717) * t631;
t572 = (t590 * t766 - t663 * t718) * t630;
t571 = (t592 * t767 + t662 * t716) * t632;
t570 = (t591 * t768 + t661 * t717) * t631;
t569 = (t590 * t769 + t660 * t718) * t630;
t568 = (t610 * t764 + (t598 * t783 - t642 * t727) * t677) * t632;
t567 = (t609 * t765 + (t596 * t784 - t642 * t730) * t676) * t631;
t566 = (t608 * t766 + (t594 * t785 - t642 * t733) * t675) * t630;
t565 = (t610 * t767 + (t597 * t783 + t642 * t728) * t677) * t632;
t564 = (t609 * t768 + (t595 * t784 + t642 * t731) * t676) * t631;
t563 = (t608 * t769 + (t593 * t785 + t642 * t734) * t675) * t630;
t1 = [m(4) + (t574 * t764 + (-t568 * t727 + t586 * t598) * t677) * t632 + (t573 * t765 + (-t567 * t730 + t584 * t596) * t676) * t631 + (t572 * t766 + (-t566 * t733 + t582 * t594) * t675) * t630, (t574 * t767 + (t568 * t728 + t586 * t597) * t677) * t632 + (t573 * t768 + (t567 * t731 + t584 * t595) * t676) * t631 + (t572 * t769 + (t566 * t734 + t582 * t593) * t675) * t630, (-t574 * t658 + t586 * t726) * t632 + (-t573 * t657 + t584 * t729) * t631 + (-t572 * t656 + t582 * t732) * t630 + (t566 * t772 + t567 * t771 + t568 * t770) * t710; (t571 * t764 + (-t565 * t727 + t585 * t598) * t677) * t632 + (t570 * t765 + (-t564 * t730 + t583 * t596) * t676) * t631 + (t569 * t766 + (-t563 * t733 + t581 * t594) * t675) * t630, m(4) + (t571 * t767 + (t565 * t728 + t585 * t597) * t677) * t632 + (t570 * t768 + (t564 * t731 + t583 * t595) * t676) * t631 + (t569 * t769 + (t563 * t734 + t581 * t593) * t675) * t630, (-t571 * t658 + t585 * t726) * t632 + (-t570 * t657 + t583 * t729) * t631 + (-t569 * t656 + t581 * t732) * t630 + (t563 * t772 + t564 * t771 + t565 * t770) * t710; (t577 * t764 + (-t580 * t727 + t589 * t598) * t677) * t632 + (t576 * t765 + (-t579 * t730 + t588 * t596) * t676) * t631 + (t575 * t766 + (-t578 * t733 + t587 * t594) * t675) * t630, (t577 * t767 + (t580 * t728 + t589 * t597) * t677) * t632 + (t576 * t768 + (t579 * t731 + t588 * t595) * t676) * t631 + (t575 * t769 + (t578 * t734 + t587 * t593) * t675) * t630, m(4) + (-t577 * t658 + t589 * t726) * t632 + (-t576 * t657 + t588 * t729) * t631 + (-t575 * t656 + t587 * t732) * t630 + (t578 * t772 + t579 * t771 + t580 * t770) * t710;];
MX  = t1;
