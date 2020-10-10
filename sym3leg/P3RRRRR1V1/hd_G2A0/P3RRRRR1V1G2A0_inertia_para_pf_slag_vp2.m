% Calculate inertia matrix for parallel robot
% P3RRRRR1V1G2A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1]';
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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 03:36
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRRRR1V1G2A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR1V1G2A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR1V1G2A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRRRR1V1G2A0_inertia_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR1V1G2A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR1V1G2A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRRRR1V1G2A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR1V1G2A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR1V1G2A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:33:55
% EndTime: 2020-08-07 03:33:57
% DurationCPUTime: 1.27s
% Computational Cost: add. (2928->274), mult. (4932->455), div. (648->8), fcn. (4014->33), ass. (0->201)
t724 = 2 * pkin(1);
t723 = 4 * Ifges(3,4);
t722 = (pkin(2) * mrSges(3,1));
t658 = qJ(2,3) + qJ(3,3);
t721 = cos(qJ(1,3) - t658) / 0.2e1 + cos(qJ(1,3) + t658) / 0.2e1;
t659 = qJ(2,2) + qJ(3,2);
t720 = cos(qJ(1,2) - t659) / 0.2e1 + cos(qJ(1,2) + t659) / 0.2e1;
t660 = qJ(2,1) + qJ(3,1);
t719 = cos(qJ(1,1) - t660) / 0.2e1 + cos(qJ(1,1) + t660) / 0.2e1;
t612 = sin(qJ(2,3));
t718 = pkin(1) * t612;
t615 = sin(qJ(2,2));
t717 = pkin(1) * t615;
t618 = sin(qJ(2,1));
t716 = pkin(1) * t618;
t715 = 2 * Ifges(2,4) - 2 * Ifges(3,4);
t607 = Ifges(3,2) - Ifges(3,1);
t611 = sin(qJ(3,3));
t714 = mrSges(3,2) * t611;
t614 = sin(qJ(3,2));
t713 = mrSges(3,2) * t614;
t617 = sin(qJ(3,1));
t712 = mrSges(3,2) * t617;
t711 = Ifges(3,4) * t611;
t710 = Ifges(3,4) * t614;
t709 = Ifges(3,4) * t617;
t632 = 0.1e1 / pkin(3);
t708 = Ifges(3,3) * t632;
t707 = m(3) * pkin(2) + mrSges(2,1);
t706 = pkin(2) * mrSges(3,3) - Ifges(2,5);
t608 = legFrame(3,2);
t594 = cos(t608);
t620 = cos(qJ(3,3));
t621 = cos(qJ(2,3));
t666 = t612 * t611;
t638 = -t620 * t621 + t666;
t663 = t621 * t611;
t639 = -t620 * t612 - t663;
t591 = sin(t608);
t613 = sin(qJ(1,3));
t678 = t591 * t613;
t549 = t639 * t594 + t638 * t678;
t601 = 0.1e1 / t611;
t705 = t549 * t601;
t675 = t594 * t613;
t550 = t639 * t591 - t638 * t675;
t704 = t550 * t601;
t609 = legFrame(2,2);
t595 = cos(t609);
t623 = cos(qJ(3,2));
t624 = cos(qJ(2,2));
t665 = t615 * t614;
t636 = -t623 * t624 + t665;
t662 = t624 * t614;
t637 = -t623 * t615 - t662;
t592 = sin(t609);
t616 = sin(qJ(1,2));
t677 = t592 * t616;
t551 = t637 * t595 + t636 * t677;
t602 = 0.1e1 / t614;
t703 = t551 * t602;
t674 = t595 * t616;
t552 = t637 * t592 - t636 * t674;
t702 = t552 * t602;
t610 = legFrame(1,2);
t596 = cos(t610);
t626 = cos(qJ(3,1));
t627 = cos(qJ(2,1));
t664 = t618 * t617;
t634 = -t626 * t627 + t664;
t661 = t627 * t617;
t635 = -t626 * t618 - t661;
t593 = sin(t610);
t619 = sin(qJ(1,1));
t676 = t593 * t619;
t553 = t635 * t596 + t634 * t676;
t603 = 0.1e1 / t617;
t701 = t553 * t603;
t673 = t596 * t619;
t554 = t635 * t593 - t634 * t673;
t700 = t554 * t603;
t588 = pkin(3) * t620 + pkin(2);
t567 = pkin(3) * t663 + t612 * t588;
t570 = -pkin(3) * t666 + t588 * t621;
t555 = t591 * t567 - t570 * t675;
t699 = t555 * t601;
t589 = pkin(3) * t623 + pkin(2);
t568 = pkin(3) * t662 + t615 * t589;
t571 = -pkin(3) * t665 + t589 * t624;
t556 = t592 * t568 - t571 * t674;
t698 = t556 * t602;
t590 = pkin(3) * t626 + pkin(2);
t569 = pkin(3) * t661 + t618 * t590;
t572 = -pkin(3) * t664 + t590 * t627;
t557 = t593 * t569 - t572 * t673;
t697 = t557 * t603;
t558 = t567 * t594 + t570 * t678;
t696 = t558 * t601;
t559 = t568 * t595 + t571 * t677;
t695 = t559 * t602;
t560 = t569 * t596 + t572 * t676;
t694 = t560 * t603;
t564 = (-Ifges(3,5) * t621 + Ifges(3,6) * t612) * t611 - t620 * (Ifges(3,5) * t612 + Ifges(3,6) * t621);
t693 = t564 * t632;
t565 = (-Ifges(3,5) * t624 + Ifges(3,6) * t615) * t614 - t623 * (Ifges(3,5) * t615 + Ifges(3,6) * t624);
t692 = t565 * t632;
t566 = (-Ifges(3,5) * t627 + Ifges(3,6) * t618) * t617 - t626 * (Ifges(3,5) * t618 + Ifges(3,6) * t627);
t691 = t566 * t632;
t622 = cos(qJ(1,3));
t690 = t570 * t622;
t625 = cos(qJ(1,2));
t689 = t571 * t625;
t628 = cos(qJ(1,1));
t688 = t572 * t628;
t576 = 0.1e1 / (pkin(3) * cos(t658) + t621 * pkin(2) + pkin(1));
t687 = t576 * t613;
t686 = t576 * t622;
t577 = 0.1e1 / (pkin(3) * cos(t659) + t624 * pkin(2) + pkin(1));
t685 = t577 * t616;
t684 = t577 * t625;
t578 = 0.1e1 / (pkin(3) * cos(t660) + t627 * pkin(2) + pkin(1));
t683 = t578 * t619;
t682 = t578 * t628;
t579 = Ifges(3,3) + (mrSges(3,1) * t620 - t714) * pkin(2);
t681 = t579 * t632;
t580 = Ifges(3,3) + (mrSges(3,1) * t623 - t713) * pkin(2);
t680 = t580 * t632;
t581 = Ifges(3,3) + (mrSges(3,1) * t626 - t712) * pkin(2);
t679 = t581 * t632;
t633 = 0.1e1 / pkin(2);
t672 = t601 * t633;
t671 = t602 * t633;
t670 = t603 * t633;
t604 = t620 ^ 2;
t669 = t607 * t604;
t605 = t623 ^ 2;
t668 = t607 * t605;
t606 = t626 ^ 2;
t667 = t607 * t606;
t599 = 2 * t722;
t598 = -0.2e1 * pkin(2) * mrSges(3,2);
t630 = m(3) * pkin(2) ^ 2;
t657 = Ifges(2,3) + Ifges(3,3) + t630;
t656 = t601 * t690;
t655 = t632 * t690;
t654 = t602 * t689;
t653 = t632 * t689;
t652 = t603 * t688;
t651 = t632 * t688;
t650 = t591 * t686;
t649 = t594 * t686;
t648 = t592 * t684;
t647 = t595 * t684;
t646 = t593 * t682;
t645 = t596 * t682;
t644 = t601 * t721;
t643 = t602 * t720;
t642 = t603 * t719;
t641 = Ifges(2,1) + Ifges(3,2) + Ifges(1,3) + (m(2) + m(3)) * (pkin(1) ^ 2);
t640 = Ifges(2,2) + t630 - Ifges(2,1) - t607;
t600 = mrSges(3,1) * t724;
t587 = t617 * t598;
t586 = t614 * t598;
t585 = t611 * t598;
t575 = t626 * t599 + t587 + t657;
t574 = t623 * t599 + t586 + t657;
t573 = t620 * t599 + t585 + t657;
t563 = (-Ifges(3,5) * t626 + t617 * Ifges(3,6) + t706) * t618 - (t617 * Ifges(3,5) + Ifges(3,6) * t626 + Ifges(2,6)) * t627;
t562 = (-Ifges(3,5) * t623 + t614 * Ifges(3,6) + t706) * t615 - (t614 * Ifges(3,5) + Ifges(3,6) * t623 + Ifges(2,6)) * t624;
t561 = (-Ifges(3,5) * t620 + t611 * Ifges(3,6) + t706) * t612 - (t611 * Ifges(3,5) + Ifges(3,6) * t620 + Ifges(2,6)) * t621;
t548 = -t566 * t683 + (-Ifges(3,3) * t651 + t581 * t719) * t670;
t547 = -t565 * t685 + (-Ifges(3,3) * t653 + t580 * t720) * t671;
t546 = -t564 * t687 + (-Ifges(3,3) * t655 + t579 * t721) * t672;
t545 = -t563 * t683 + (t575 * t719 - t581 * t651) * t670;
t544 = -t562 * t685 + (t574 * t720 - t580 * t653) * t671;
t543 = -t561 * t687 + (t573 * t721 - t579 * t655) * t672;
t542 = -t667 + 0.2e1 * (-mrSges(3,2) * t716 - t709) * t626 - 0.2e1 * (t617 * mrSges(3,1) + mrSges(2,2)) * t716 + t641 + (t600 * t626 + (t707 - t712) * t724 + (t606 * t723 + t598 * t626 + 0.2e1 * (-t607 * t626 - t722) * t617 + t715) * t618 + (0.2e1 * t667 + (t599 + 0.4e1 * t709) * t626 + t587 + t640) * t627) * t627;
t541 = -t668 + 0.2e1 * (-mrSges(3,2) * t717 - t710) * t623 - 0.2e1 * (t614 * mrSges(3,1) + mrSges(2,2)) * t717 + t641 + (t600 * t623 + (t707 - t713) * t724 + (t605 * t723 + t598 * t623 + 0.2e1 * (-t607 * t623 - t722) * t614 + t715) * t615 + (0.2e1 * t668 + (t599 + 0.4e1 * t710) * t623 + t586 + t640) * t624) * t624;
t540 = -t669 + 0.2e1 * (-mrSges(3,2) * t718 - t711) * t620 - 0.2e1 * (t611 * mrSges(3,1) + mrSges(2,2)) * t718 + t641 + (t600 * t620 + (t707 - t714) * t724 + (t604 * t723 + t598 * t620 + 0.2e1 * (-t607 * t620 - t722) * t611 + t715) * t612 + (0.2e1 * t669 + (t599 + 0.4e1 * t711) * t620 + t585 + t640) * t621) * t621;
t539 = t566 * t645 + (t554 * t581 + t557 * t708) * t670;
t538 = t565 * t647 + (t552 * t580 + t556 * t708) * t671;
t537 = t564 * t649 + (t550 * t579 + t555 * t708) * t672;
t536 = -t566 * t646 + (t553 * t581 + t560 * t708) * t670;
t535 = -t565 * t648 + (t551 * t580 + t559 * t708) * t671;
t534 = -t564 * t650 + (t549 * t579 + t558 * t708) * t672;
t533 = t563 * t645 + (t554 * t575 + t557 * t679) * t670;
t532 = t562 * t647 + (t552 * t574 + t556 * t680) * t671;
t531 = t561 * t649 + (t550 * t573 + t555 * t681) * t672;
t530 = -t563 * t646 + (t553 * t575 + t560 * t679) * t670;
t529 = -t562 * t648 + (t551 * t574 + t559 * t680) * t671;
t528 = -t561 * t650 + (t549 * t573 + t558 * t681) * t672;
t527 = -t542 * t683 + (t563 * t719 - t566 * t651) * t670;
t526 = -t541 * t685 + (t562 * t720 - t565 * t653) * t671;
t525 = -t540 * t687 + (t561 * t721 - t564 * t655) * t672;
t524 = t542 * t645 + (t554 * t563 + t557 * t691) * t670;
t523 = t541 * t647 + (t552 * t562 + t556 * t692) * t671;
t522 = t540 * t649 + (t550 * t561 + t555 * t693) * t672;
t521 = -t542 * t646 + (t553 * t563 + t560 * t691) * t670;
t520 = -t541 * t648 + (t551 * t562 + t559 * t692) * t671;
t519 = -t540 * t650 + (t549 * t561 + t558 * t693) * t672;
t1 = [t522 * t649 + t523 * t647 + t524 * t645 + m(4) + (t531 * t704 + t532 * t702 + t533 * t700 + (t537 * t699 + t538 * t698 + t539 * t697) * t632) * t633, -t522 * t650 - t523 * t648 - t524 * t646 + (t531 * t705 + t532 * t703 + t533 * t701 + (t537 * t696 + t538 * t695 + t539 * t694) * t632) * t633, -t522 * t687 - t523 * t685 - t524 * t683 + (t533 * t642 + t532 * t643 + t531 * t644 + (-t537 * t656 - t538 * t654 - t539 * t652) * t632) * t633; t519 * t649 + t520 * t647 + t521 * t645 + (t528 * t704 + t529 * t702 + t530 * t700 + (t534 * t699 + t535 * t698 + t536 * t697) * t632) * t633, -t519 * t650 - t520 * t648 - t521 * t646 + m(4) + (t528 * t705 + t529 * t703 + t530 * t701 + (t534 * t696 + t535 * t695 + t536 * t694) * t632) * t633, -t519 * t687 - t520 * t685 - t521 * t683 + (t530 * t642 + t529 * t643 + t528 * t644 + (-t534 * t656 - t535 * t654 - t536 * t652) * t632) * t633; t525 * t649 + t526 * t647 + t527 * t645 + (t543 * t704 + t544 * t702 + t545 * t700 + (t546 * t699 + t547 * t698 + t548 * t697) * t632) * t633, -t525 * t650 - t526 * t648 - t527 * t646 + (t543 * t705 + t544 * t703 + t545 * t701 + (t546 * t696 + t547 * t695 + t548 * t694) * t632) * t633, -t525 * t687 - t526 * t685 - t527 * t683 + m(4) + (t545 * t642 + t544 * t643 + t543 * t644 + (-t546 * t656 - t547 * t654 - t548 * t652) * t632) * t633;];
MX  = t1;
