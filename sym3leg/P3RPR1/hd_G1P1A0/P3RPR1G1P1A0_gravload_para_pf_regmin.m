% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPR1G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x10]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:58
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPR1G1P1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPR1G1P1A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RPR1G1P1A0_gravload_para_pf_regmin: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPR1G1P1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RPR1G1P1A0_gravload_para_pf_regmin: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPR1G1P1A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPR1G1P1A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:58:13
% EndTime: 2019-05-03 14:58:14
% DurationCPUTime: 0.56s
% Computational Cost: add. (479->116), mult. (862->220), div. (54->3), fcn. (748->14), ass. (0->102)
t628 = legFrame(3,3);
t620 = sin(t628);
t623 = cos(t628);
t588 = g(1) * t620 - g(2) * t623;
t591 = g(1) * t623 + g(2) * t620;
t631 = sin(qJ(1,3));
t634 = cos(qJ(1,3));
t564 = t588 * t631 - t591 * t634;
t629 = legFrame(2,3);
t621 = sin(t629);
t624 = cos(t629);
t589 = g(1) * t621 - g(2) * t624;
t592 = g(1) * t624 + g(2) * t621;
t632 = sin(qJ(1,2));
t635 = cos(qJ(1,2));
t566 = t589 * t632 - t592 * t635;
t630 = legFrame(1,3);
t622 = sin(t630);
t625 = cos(t630);
t590 = g(1) * t622 - g(2) * t625;
t593 = g(1) * t625 + g(2) * t622;
t633 = sin(qJ(1,1));
t636 = cos(qJ(1,1));
t568 = t590 * t633 - t593 * t636;
t587 = -t633 * t622 + t636 * t625;
t642 = 0.1e1 / qJ(2,1);
t649 = t587 * t642;
t585 = -t632 * t621 + t635 * t624;
t641 = 0.1e1 / qJ(2,2);
t651 = t585 * t641;
t583 = -t631 * t620 + t634 * t623;
t640 = 0.1e1 / qJ(2,3);
t653 = t583 * t640;
t661 = t564 * t653 + t566 * t651 + t568 * t649;
t586 = t636 * t622 + t633 * t625;
t650 = t586 * t642;
t584 = t635 * t621 + t632 * t624;
t652 = t584 * t641;
t582 = t634 * t620 + t631 * t623;
t654 = t582 * t640;
t660 = t564 * t654 + t566 * t652 + t568 * t650;
t645 = koppelP(1,2);
t648 = koppelP(1,1);
t600 = t633 * t648 - t636 * t645;
t601 = t633 * t645 + t636 * t648;
t639 = xP(3);
t626 = sin(t639);
t627 = cos(t639);
t560 = (t600 * t627 - t601 * t626) * t625 + t622 * (t600 * t626 + t601 * t627);
t655 = t560 * t642;
t644 = koppelP(2,2);
t647 = koppelP(2,1);
t598 = t632 * t647 - t635 * t644;
t599 = t632 * t644 + t635 * t647;
t559 = (t598 * t627 - t599 * t626) * t624 + t621 * (t598 * t626 + t599 * t627);
t656 = t559 * t641;
t643 = koppelP(3,2);
t646 = koppelP(3,1);
t596 = t631 * t646 - t634 * t643;
t597 = t631 * t643 + t634 * t646;
t558 = (t596 * t627 - t597 * t626) * t623 + t620 * (t596 * t626 + t597 * t627);
t657 = t558 * t640;
t659 = t564 * t657 + t566 * t656 + t568 * t655;
t658 = pkin(1) * g(2);
t565 = t588 * t634 + t591 * t631;
t567 = t589 * t635 + t592 * t632;
t569 = t590 * t636 + t593 * t633;
t638 = pkin(1) * g(1);
t637 = pkin(1) + pkin(2);
t619 = g(1) * qJ(2,1) + t658;
t618 = -g(2) * qJ(2,1) + t638;
t617 = g(1) * qJ(2,2) + t658;
t616 = -g(2) * qJ(2,2) + t638;
t615 = g(1) * qJ(2,3) + t658;
t614 = -g(2) * qJ(2,3) + t638;
t613 = qJ(2,1) * t648 + t637 * t645;
t612 = qJ(2,2) * t647 + t637 * t644;
t611 = qJ(2,3) * t646 + t637 * t643;
t610 = -qJ(2,1) * t645 + t637 * t648;
t609 = -qJ(2,2) * t644 + t637 * t647;
t608 = -qJ(2,3) * t643 + t637 * t646;
t607 = qJ(2,1) * t633 + t636 * t637;
t606 = qJ(2,2) * t632 + t635 * t637;
t605 = qJ(2,3) * t631 + t634 * t637;
t604 = -qJ(2,1) * t636 + t633 * t637;
t603 = -qJ(2,2) * t635 + t632 * t637;
t602 = -qJ(2,3) * t634 + t631 * t637;
t595 = g(1) * t627 + g(2) * t626;
t594 = g(1) * t626 - g(2) * t627;
t581 = t610 * t633 - t613 * t636;
t580 = t609 * t632 - t612 * t635;
t579 = t608 * t631 - t611 * t634;
t578 = t610 * t636 + t613 * t633;
t577 = t609 * t635 + t612 * t632;
t576 = t608 * t634 + t611 * t631;
t563 = (t618 * t633 - t619 * t636) * t625 + (t618 * t636 + t619 * t633) * t622;
t562 = (t616 * t632 - t617 * t635) * t624 + (t616 * t635 + t617 * t632) * t621;
t561 = (t614 * t631 - t615 * t634) * t623 + (t614 * t634 + t615 * t631) * t620;
t557 = t565 * t653 + t567 * t651 + t569 * t649;
t556 = t565 * t654 + t567 * t652 + t569 * t650;
t555 = t565 * t657 + t567 * t656 + t569 * t655;
t1 = [0, t557, -t661, t557, t661, (t587 * t563 - (-t604 * t622 + t607 * t625) * t569) * t642 + (t585 * t562 - (-t603 * t621 + t606 * t624) * t567) * t641 + (t583 * t561 - (-t602 * t620 + t605 * t623) * t565) * t640, 0, 0, 0, -t594 * t626 - t595 * t627; 0, t556, -t660, t556, t660, (t586 * t563 - (t604 * t625 + t607 * t622) * t569) * t642 + (t584 * t562 - (t603 * t624 + t606 * t621) * t567) * t641 + (t582 * t561 - (t602 * t623 + t605 * t620) * t565) * t640, 0, 0, 0, t594 * t627 - t595 * t626; 0, t555, -t659, t555, t659, (t560 * t563 - ((-t578 * t626 + t581 * t627) * t625 + (t578 * t627 + t581 * t626) * t622) * t569) * t642 + (t559 * t562 - ((-t577 * t626 + t580 * t627) * t624 + (t577 * t627 + t580 * t626) * t621) * t567) * t641 + (t558 * t561 - ((-t576 * t626 + t579 * t627) * t623 + (t576 * t627 + t579 * t626) * t620) * t565) * t640, 0, t594, t595, 0;];
tau_reg  = t1;
