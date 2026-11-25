# 전자공학 시각화 튜터

전자공학 전공 공부용 시각 자료 모음입니다. 각 파일은 서버 없이 단독으로 동작하는 HTML 한 개입니다.

## 도구

| 파일 | 과목 | 내용 |
|---|---|---|
| `convolution.html` | 신호및시스템 | 컨볼루션 4단계(뒤집기·이동·곱하기·적분), τ 슬라이더, 해석해 대조 |
| `eigen.html` | 선형대수 | 2×2 고유값·고유벡터, 단위원→타원 변환, 고유방향 정렬 감지, 복소 고유값 |
| `gauss.html` | 전자장 | 가우스 법칙 — 점·선·평면·구 대칭, 가우스면 조작, 풀이 과정 |
| `pipeline.html` | 컴퓨터구조 | 5단 파이프라인 RAW 해저드, 스톨 버블, 포워딩 ON/OFF, CPI |
| `mosfet.html` | 반도체소자 | 컷오프/선형/포화 판정, 출력·전달 특성, 채널 에너지(핀치오프), g_m·r_o |
| `mos-cap.html` | 반도체소자 | 축적/평탄/공핍/반전 판정, MOS C-V, 게이트–산화막 개략 |
| `schottky-band.html` | 반도체소자 | Schottky 밴드 다이어그램 (평형/역바이어스) |
| `pmos-cmos.html` | 반도체소자 | p채널 MOSFET · CMOS 인버터 VTC와 스위칭 문턱 |
| `three-phase.html` | 회로이론 | Y-Y·Y-Δ·Δ-Y·Δ-Δ, 평형/불평형, 4선식/3선식, 페이저도, 풀이 과정 |
| `circuit-editor.html` | 회로이론 | 회로를 직접 그리고 테브난/노턴 등가 구하기 (교류 복소 임피던스 포함) |
| `thevenin.html` | 회로이론 | 값 편집형 테브난/노턴 단계 확인 |

`index.html` 이 목록 페이지입니다.

## 계산 검증

각 도구는 Python 스크립트로 독립 계산한 값과 대조했습니다.

- `verify_mosfet.py` — 4개 동작점, 경계 연속성, g_m·r_o 수치미분
- `verify_mosc.py` — MOS-C 7개 항목
- `verify_tyu91.py` — 교재 Schottky 예제 6개 값 + 불변식 5개
- `verify_3phase.py` — 4가지 결선 기준값
- `verify_thevenin.py` — 세 가지 경로로 구한 R_th 대조
- `verify_ac_thevenin.py` — 코일·콘덴서 교류 Z_th = V_th/I_sc 대조
- `verify_pmos_cmos.py` — pMOS·CMOS 인버터 기준값 대조
- `verify_convolution.py` — 이산 합 vs 해석해(사다리꼴·구간별 닫힌형), 교환법칙, δ 항등식, 격자 1차 수렴
- `verify_eigen.py` — 특성다항식 잔차, ‖Av−λv‖, Vieta, 거듭제곱법, 복소 회전각
- `verify_gauss.py` — 구면 플럭스 직접 적분 = Q/ε₀, 껍질 상쇄, 몬테카를로 Q_enc, 선·평면 쿨롱 합 대조
- `verify_pipeline.py` — 손계산 스케줄 8케이스 + 구조 항등식, 해저드 원장부($0·sw·최근 생산자)

실행: `python verify_mosfet.py` (종료 코드 0이면 통과)

## 범위와 한계

- MOSFET은 1차 square-law 모델입니다 (subthreshold, 속도포화, DIBL 없음).
- 3상은 전원이 평형일 때만 다룹니다. 부하는 상별 불평형까지 지원합니다.
- 테브난 도구는 독립전원만, 직병렬 축소만 지원합니다.
- 컨볼루션은 격자 합(Δ) 기반이라 신호 경계에서 O(Δ) 오차가 있습니다 (페이지 내 검산 참조).
- 파이프라인은 RAW 데이터 해저드만 다룹니다. 레지스터는 WB 다음 사이큼에 읽는 관례를 쓰며, 분기·구조 해저드는 미지원입니다.
- 가우스 도구의 선·평면 전하는 유한 크기 근사로 검증했습니다 (무한 극한과 대조).

개인 학습용으로 만든 자료입니다.
