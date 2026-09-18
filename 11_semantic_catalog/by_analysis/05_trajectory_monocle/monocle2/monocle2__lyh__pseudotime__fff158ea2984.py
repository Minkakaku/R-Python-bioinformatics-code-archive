import streamlit as st
from utils.utils import *


def pseudotime():
    _, context, _ = st.columns([0.5, 8, 1])
    _.button(label='↩ 返回', on_click=callback, type='primary')
    st.write('pseudotime')

