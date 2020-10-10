% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRR1G1A0
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x8]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRR1G1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G1A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G1A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR1G1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G1A0_gravload_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G1A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G1A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:14:56
% EndTime: 2020-03-09 21:14:56
% DurationCPUTime: 0.26s
% Computational Cost: add. (739->76), mult. (457->139), div. (84->5), fcn. (552->18), ass. (0->75)
t622 = legFrame(3,3);
t613 = sin(t622);
t616 = cos(t622);
t589 = g(1) * t613 - t616 * g(2);
t592 = g(1) * t616 + g(2) * t613;
t619 = pkin(7) + qJ(2,3);
t610 = qJ(3,3) + t619;
t595 = sin(t610);
t598 = cos(t610);
t568 = t589 * t598 + t592 * t595;
t604 = sin(t619);
t607 = cos(t619);
t580 = 0.1e1 / (-t595 * t607 + t598 * t604);
t641 = t568 * t580;
t623 = legFrame(2,3);
t614 = sin(t623);
t617 = cos(t623);
t590 = g(1) * t614 - t617 * g(2);
t593 = g(1) * t617 + g(2) * t614;
t620 = pkin(7) + qJ(2,2);
t611 = qJ(3,2) + t620;
t596 = sin(t611);
t599 = cos(t611);
t569 = t590 * t599 + t593 * t596;
t605 = sin(t620);
t608 = cos(t620);
t581 = 0.1e1 / (-t596 * t608 + t599 * t605);
t640 = t569 * t581;
t624 = legFrame(1,3);
t615 = sin(t624);
t618 = cos(t624);
t591 = g(1) * t615 - t618 * g(2);
t594 = g(1) * t618 + g(2) * t615;
t621 = pkin(7) + qJ(2,1);
t612 = qJ(3,1) + t621;
t597 = sin(t612);
t600 = cos(t612);
t570 = t591 * t600 + t594 * t597;
t606 = sin(t621);
t609 = cos(t621);
t582 = 0.1e1 / (-t597 * t609 + t600 * t606);
t639 = t570 * t582;
t571 = -t589 * t595 + t592 * t598;
t638 = t571 * t580;
t572 = -t590 * t596 + t593 * t599;
t637 = t572 * t581;
t573 = -t591 * t597 + t594 * t600;
t636 = t573 * t582;
t629 = t595 * t616 + t598 * t613;
t635 = t580 * t629;
t584 = t595 * t613 - t598 * t616;
t634 = t580 * t584;
t628 = t596 * t617 + t599 * t614;
t633 = t581 * t628;
t586 = t596 * t614 - t599 * t617;
t632 = t581 * t586;
t627 = t597 * t618 + t600 * t615;
t631 = t582 * t627;
t588 = t597 * t615 - t600 * t618;
t630 = t582 * t588;
t626 = 0.1e1 / pkin(2);
t625 = 0.1e1 / pkin(3);
t579 = -t591 * t606 + t594 * t609;
t578 = -t590 * t605 + t593 * t608;
t577 = -t589 * t604 + t592 * t607;
t576 = t591 * t609 + t594 * t606;
t575 = t590 * t608 + t593 * t605;
t574 = t589 * t607 + t592 * t604;
t567 = -pkin(2) * (t606 * t615 - t609 * t618) - t588 * pkin(3);
t566 = -pkin(2) * (t605 * t614 - t608 * t617) - t586 * pkin(3);
t565 = -pkin(2) * (t604 * t613 - t607 * t616) - t584 * pkin(3);
t564 = pkin(2) * (t606 * t618 + t609 * t615) + t627 * pkin(3);
t563 = pkin(2) * (t605 * t617 + t608 * t614) + t628 * pkin(3);
t562 = pkin(2) * (t604 * t616 + t613 * t607) + t629 * pkin(3);
t1 = [0, 0, (t574 * t634 + t575 * t632 + t576 * t630) * t626, (t577 * t634 + t578 * t632 + t579 * t630) * t626, 0, (t568 * t634 + t569 * t632 + t570 * t630 + (t565 * t641 + t566 * t640 + t567 * t639) * t625) * t626, (t571 * t634 + t572 * t632 + t573 * t630 + (t565 * t638 + t566 * t637 + t567 * t636) * t625) * t626, -g(1); 0, 0, (-t574 * t635 - t575 * t633 - t576 * t631) * t626, (-t577 * t635 - t578 * t633 - t579 * t631) * t626, 0, (-t568 * t635 - t569 * t633 - t570 * t631 + (t562 * t641 + t563 * t640 + t564 * t639) * t625) * t626, (-t571 * t635 - t572 * t633 - t573 * t631 + (t562 * t638 + t563 * t637 + t564 * t636) * t625) * t626, -g(2); -3 * g(3), 0, 0, 0, 0, 0, 0, -g(3);];
tau_reg  = t1;
